export type ApiErrorPayload = {
  error?: string
}

const flattenErrorValue = (value: unknown): string[] => {
  if (typeof value === 'string') {
    return [value]
  }

  if (Array.isArray(value)) {
    return value.flatMap((entry) => flattenErrorValue(entry))
  }

  if (!value || typeof value !== 'object') {
    return []
  }

  return Object.entries(value as Record<string, unknown>).flatMap(([key, entryValue]) => {
    const nestedMessages = flattenErrorValue(entryValue)

    if (nestedMessages.length === 0) {
      return []
    }

    return nestedMessages.map((message) => `${key}: ${message}`)
  })
}

const getDataError = (obj: Record<string, unknown>): string | null => {
  if (!('data' in obj) || !obj.data) return null

  const data = obj.data as ApiErrorPayload
  return typeof data.error === 'string' ? data.error : null
}

const getDataValidationError = (obj: Record<string, unknown>): string | null => {
  if (!('data' in obj) || !obj.data || typeof obj.data !== 'object') return null

  const messages = flattenErrorValue(obj.data)
  return messages.length > 0 ? messages.join(' | ') : null
}

const getResponseError = (obj: Record<string, unknown>): string | null => {
  if (!('response' in obj) || !obj.response || typeof obj.response !== 'object') return null

  const response = obj.response as Record<string, unknown>
  if (!('_data' in response) || !response._data) return null

  const data = response._data as ApiErrorPayload
  return typeof data.error === 'string' ? data.error : null
}

const getResponseValidationError = (obj: Record<string, unknown>): string | null => {
  if (!('response' in obj) || !obj.response || typeof obj.response !== 'object') return null

  const response = obj.response as Record<string, unknown>
  if (!('_data' in response) || !response._data || typeof response._data !== 'object') return null

  const messages = flattenErrorValue(response._data)
  return messages.length > 0 ? messages.join(' | ') : null
}

const getMessageError = (obj: Record<string, unknown>): string | null => {
  return 'message' in obj && typeof obj.message === 'string' ? obj.message : null
}

export const getErrorMessage = (error: unknown): string => {
  if (!error || typeof error !== 'object') return String(error)

  const obj = error as Record<string, unknown>
  return (
    getDataError(obj) ??
    getResponseError(obj) ??
    getDataValidationError(obj) ??
    getResponseValidationError(obj) ??
    getMessageError(obj) ??
    String(error)
  )
}
