export type ApiErrorPayload = {
  error?: string
}

const getDataError = (obj: Record<string, unknown>): string | null => {
  if (!('data' in obj) || !obj.data) return null

  const data = obj.data as ApiErrorPayload
  return typeof data.error === 'string' ? data.error : null
}

const getResponseError = (obj: Record<string, unknown>): string | null => {
  if (!('response' in obj) || !obj.response || typeof obj.response !== 'object') return null

  const response = obj.response as Record<string, unknown>
  if (!('_data' in response) || !response._data) return null

  const data = response._data as ApiErrorPayload
  return typeof data.error === 'string' ? data.error : null
}

const getMessageError = (obj: Record<string, unknown>): string | null => {
  return 'message' in obj && typeof obj.message === 'string' ? obj.message : null
}

export const getErrorMessage = (error: unknown): string => {
  if (!error || typeof error !== 'object') return String(error)

  const obj = error as Record<string, unknown>
  return getDataError(obj) ?? getResponseError(obj) ?? getMessageError(obj) ?? String(error)
}
