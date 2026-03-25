import { defineStore } from 'pinia'

export type StreamEvent =
  | { type: 'content_delta'; delta: string }
  | { type: 'thinking_delta'; delta: string }
  | { type: 'done'; duration_ms: number }
  | { type: 'error'; error: string }

export type StreamHandlers = {
  onContentDelta: (delta: string) => void
  onThinkingDelta: (delta: string) => void
  onDone: (durationMs: number) => void
  onError: (error: string) => void
}

const parseSseLine = (line: string): StreamEvent | null => {
  if (!line.startsWith('data:')) return null

  const payloadText = line.slice(5).trim()
  if (!payloadText) return null

  try {
    return JSON.parse(payloadText) as StreamEvent
  } catch {
    return null
  }
}

const SSE_EVENT_SEPARATOR_REGEX = /\r?\n\r?\n/

export const useMlStreamingStore = defineStore('mlStreamingStore', {
  state: () => ({}),

  actions: {
    async processStream(reader: ReadableStreamDefaultReader<Uint8Array>, handlers: StreamHandlers) {
      const decoder = new TextDecoder()
      let buffer = ''

      for (;;) {
        const { value, done } = await reader.read()
        if (done) break

        buffer += decoder.decode(value, { stream: true })
        buffer = this._drainSseBuffer(buffer, handlers)
      }

      this._drainSseBuffer(buffer, handlers)
    },

    _drainSseBuffer(buffer: string, handlers: StreamHandlers): string {
      const MAX_DRAIN_CYCLES = 100

      for (let i = 0; i < MAX_DRAIN_CYCLES; i++) {
        const nextBuffer = this._processSseBuffer(buffer, handlers)
        if (nextBuffer === buffer) return buffer
        buffer = nextBuffer
      }

      return buffer
    },

    _processSseBuffer(buffer: string, handlers: StreamHandlers): string {
      for (let i = 0; i < 100; i++) {
        const separatorMatch = buffer.match(SSE_EVENT_SEPARATOR_REGEX)
        if (!separatorMatch || separatorMatch.index == null) break

        const chunk = buffer.slice(0, separatorMatch.index)
        buffer = buffer.slice(separatorMatch.index + separatorMatch[0].length)

        this._processSseChunk(chunk, handlers)
      }

      return buffer
    },

    _processSseChunk(chunk: string, handlers: StreamHandlers) {
      const lines = chunk.split('\n')
      for (const line of lines) {
        const payload = parseSseLine(line)
        if (payload) {
          this._processStreamEvent(payload, handlers)
        }
      }
    },

    _processStreamEvent(event: StreamEvent, handlers: StreamHandlers) {
      switch (event.type) {
        case 'content_delta':
          if (event.delta) handlers.onContentDelta(event.delta)
          break
        case 'thinking_delta':
          if (event.delta) handlers.onThinkingDelta(event.delta)
          break
        case 'done':
          handlers.onDone(event.duration_ms)
          break
        case 'error':
          handlers.onError(event.error)
          break
      }
    },
  },
})
