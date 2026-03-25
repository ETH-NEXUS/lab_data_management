import { defineStore } from 'pinia'
import type { ChatMessage, ChatStatus } from '~~/types/ml'
import { getErrorMessage } from 'errors'

import { useMlModelsStore } from './models'
import { useMlStreamingStore } from './streaming'

const isAbortError = (error: unknown): boolean => {
  if (!error || typeof error !== 'object') return false
  return 'name' in error && (error as { name?: unknown }).name === 'AbortError'
}

const getResponseErrorText = async (response: Response): Promise<string | null> => {
  const contentType = response.headers.get('content-type') || ''

  try {
    if (contentType.includes('application/json')) {
      const data = (await response.json()) as unknown
      if (data && typeof data === 'object') {
        const obj = data as Record<string, unknown>
        const message = obj.error ?? obj.detail ?? obj.message
        if (typeof message === 'string' && message.trim()) return message.trim()
      }
      return JSON.stringify(data)
    }

    const text = await response.text()
    return text.trim() || null
  } catch {
    return null
  }
}

const getTextFromParts = (parts: ChatMessage['parts']): string => {
  return parts
    .filter((part) => part.type === 'text')
    .map((part) => part.text)
    .join('')
}

export interface MlChatState {
  input: string
  status: ChatStatus
  messages: ChatMessage[]
  requestStartedAtMs: number | null
  lastDurationMs: number | null
  currentTimeMs: number
  _tick: ReturnType<typeof setInterval> | null
  _activeRequestId: number
  _abortController: AbortController | null
}

export const useMlChatStore = defineStore('mlChatStore', {
  state: (): MlChatState => ({
    input: '',
    status: 'ready' as ChatStatus,
    messages: [],
    requestStartedAtMs: null as number | null,
    lastDurationMs: null as number | null,
    currentTimeMs: 0,
    _tick: null as ReturnType<typeof setInterval> | null,
    _activeRequestId: 0,
    _abortController: null as AbortController | null,
  }),

  getters: {
    durationLabel(state): string {
      const ms =
        state.requestStartedAtMs == null
          ? state.lastDurationMs
          : Math.max(0, Math.round(state.currentTimeMs - state.requestStartedAtMs))

      if (ms == null) return ''

      const s = ms / 1000
      return s < 10 ? `${s.toFixed(2)}s` : `${s.toFixed(1)}s`
    },
  },

  actions: {
    abortActiveRequest() {
      this._activeRequestId += 1

      if (this._abortController) {
        this._abortController.abort()
        this._abortController = null
      }

      if (this.status === 'submitted' || this.status === 'streaming') {
        this.status = 'ready'
      }

      this.requestStartedAtMs = null
      this._stopTimer()
    },

    _startTimer() {
      if (this._tick) return

      this.currentTimeMs = performance.now()

      this._tick = setInterval(() => {
        this.currentTimeMs = performance.now()
      }, 100)
    },

    _stopTimer() {
      if (this._tick) clearInterval(this._tick)
      this._tick = null
    },

    _createUserMessage(content: string): ChatMessage {
      const modelsStore = useMlModelsStore()
      return {
        id: `msg-${Date.now()}`,
        role: 'user',
        parts: [{ type: 'text', text: content }],
        metadata: { model: modelsStore.selectedModel },
      }
    },

    _createAssistantMessage(): ChatMessage {
      const modelsStore = useMlModelsStore()
      return {
        id: `msg-${Date.now()}-assistant`,
        role: 'assistant',
        parts: [{ type: 'text', text: '' }],
        metadata: { model: modelsStore.selectedModel, durationMs: null, thinking: '' },
      }
    },

    _createErrorMessage(error: string, durationMs: number | null): ChatMessage {
      return {
        id: `msg-${Date.now()}-error`,
        role: 'assistant',
        parts: [{ type: 'text', text: `Error: ${error}` }],
        metadata: { error, durationMs },
      }
    },

    _getPayloadMessages() {
      return this.messages
        .filter((message) => !message.metadata?.error)
        .map((message) => ({ role: message.role, content: getTextFromParts(message.parts).trim() }))
        .filter((message) => message.content.length > 0)
    },

    _ensureAssistantMessage(assistantMessage: ChatMessage) {
      const isAdded = this.messages.some((message) => message.id === assistantMessage.id)
      if (!isAdded) {
        this.messages = [...this.messages, assistantMessage]
      }
    },

    _triggerReactivity() {
      this.messages = [...this.messages]
    },

    _handleContentDelta(assistantMessage: ChatMessage, delta: string) {
      this._ensureAssistantMessage(assistantMessage)
      assistantMessage.parts[0]!.text += delta
      this._triggerReactivity()
    },

    _handleThinkingDelta(assistantMessage: ChatMessage, delta: string) {
      this._ensureAssistantMessage(assistantMessage)

      const existing = assistantMessage.metadata?.thinking || ''

      assistantMessage.metadata = {
        ...assistantMessage.metadata,
        thinking: existing + delta,
      }

      this._triggerReactivity()
    },

    _handleStreamDone(assistantMessage: ChatMessage, durationMs: number) {
      this._ensureAssistantMessage(assistantMessage)
      this.lastDurationMs = durationMs

      assistantMessage.metadata = {
        ...assistantMessage.metadata,
        durationMs,
      }

      const i = this.messages.findIndex((message) => message.id === assistantMessage.id)
      if (i !== -1) {
        this.messages[i] = { ...assistantMessage }
      }

      this.status = 'ready'
    },

    _handleStreamError(assistantMessage: ChatMessage, error: string) {
      this._ensureAssistantMessage(assistantMessage)
      this.status = 'error'

      assistantMessage.parts[0]!.text = `Error: ${error}`

      assistantMessage.metadata = {
        ...assistantMessage.metadata,
        error,
      }
    },

    async submit() {
      const content = this.input.trim()
      if (!content) return

      this.abortActiveRequest()
      const requestId = this._activeRequestId

      this._initializeSubmission(content)

      try {
        await this._executeSubmission(requestId)
      } catch (error) {
        if (this._activeRequestId !== requestId) return
        if (isAbortError(error)) return
        this._handleSubmissionError(error)
      } finally {
        if (this._activeRequestId === requestId) {
          this._cleanupSubmission()
        }
      }
    },

    _initializeSubmission(content: string) {
      this.status = 'submitted'
      this.requestStartedAtMs = performance.now()
      this._startTimer()
      this.messages = [...this.messages, this._createUserMessage(content)]
      this.input = ''
    },

    async _executeSubmission(requestId: number) {
      const modelsStore = useMlModelsStore()
      if (!modelsStore.selectedModel) {
        throw new Error('Please select a model')
      }

      this._abortController = new AbortController()

      const response = await fetch('/api/v1/ml/chat/', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json', Accept: 'text/event-stream' },
        body: JSON.stringify({
          model: modelsStore.selectedModel,
          messages: this._getPayloadMessages(),
        }),
        signal: this._abortController.signal,
      })

      if (this._activeRequestId !== requestId) return

      if (!response.ok) {
        const details = await getResponseErrorText(response)
        const message = details || `Request failed (${response.status} ${response.statusText})`
        throw new Error(message)
      }

      const contentType = response.headers.get('content-type') || ''
      if (!contentType.includes('text/event-stream')) {
        const details = await getResponseErrorText(response)
        const suffix = details ? `: ${details}` : ''
        throw new Error(`Expected text/event-stream but got ${contentType || '(missing content-type)'}${suffix}`)
      }

      if (!response.body) {
        throw new Error('Streaming response not available')
      }

      this.status = 'streaming'
      const reader = response.body.getReader()
      const assistantMessage = this._createAssistantMessage()

      const streamingStore = useMlStreamingStore()
      await streamingStore.processStream(reader, {
        onContentDelta: (delta) => {
          if (this._activeRequestId !== requestId) return
          this._handleContentDelta(assistantMessage, delta)
        },
        onThinkingDelta: (delta) => {
          if (this._activeRequestId !== requestId) return
          this._handleThinkingDelta(assistantMessage, delta)
        },
        onDone: (durationMs) => {
          if (this._activeRequestId !== requestId) return
          this._handleStreamDone(assistantMessage, durationMs)
        },
        onError: (error) => {
          if (this._activeRequestId !== requestId) return
          this._handleStreamError(assistantMessage, error)
        },
      })
    },

    _handleSubmissionError(error: unknown) {
      this.lastDurationMs =
        this.requestStartedAtMs == null ? null : Math.max(0, Math.round(performance.now() - this.requestStartedAtMs))
      this.status = 'error'
      const errorMessage = getErrorMessage(error)
      this.messages = [...this.messages, this._createErrorMessage(errorMessage, this.lastDurationMs)]
    },

    _cleanupSubmission() {
      this.requestStartedAtMs = null
      this._abortController = null
      this._stopTimer()
    },
  },
})
