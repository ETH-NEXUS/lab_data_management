export type ChatStatus = 'submitted' | 'streaming' | 'ready' | 'error'

export type ChatMessageMetadata = {
  model?: string | null
  durationMs?: number | null
  thinking?: string
  error?: string
}

export type ChatMessage = {
  id: string
  role: 'user' | 'assistant' | 'system'
  parts: [{ type: 'text'; text: string }]
  metadata?: ChatMessageMetadata
}

export type ChatResponse = {
  content: string
  duration_ms: number
  thinking?: string
  error?: string
}

export type ModelItem = { label: string; value: string | null; icon?: string }
