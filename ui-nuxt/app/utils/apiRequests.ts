import { useAPI } from '~/composables/useAPI'

type ApiRequestMethod = 'GET' | 'POST' | 'PUT' | 'DELETE' | 'PATCH' | 'HEAD' | 'OPTIONS'

type ApiRequestOptions = {
  method?: ApiRequestMethod
  headers?: Record<string, string>
  body?: unknown
  params?: Record<string, string>
  responseType?: 'json' | 'text' | 'blob' | 'arrayBuffer'
}

/**
 * Requests one endpoint and returns response body data.
 *
 * Accepted input examples:
 * - `endpoint = 'projects/'`, options = `{ method: 'POST', body: { name: 'Screening 2026' } }`
 * - `endpoint = 'experiments/21/'`, options = `{ method: 'PATCH', body: { name: 'Dose response v2' } }`
 *
 * Returned data examples:
 * - `{ id: 21, name: 'Dose response v2' }`
 * - `{ id: 44, name: 'Screening 2026' }`
 */
export const requestApiData = async <T>(
  endpoint: string,
  options: ApiRequestOptions,
  fallbackErrorMessage: string,
): Promise<T> => {
  const { data, error } = await useAPI<T>(endpoint, options)

  if (error.value || !data.value) {
    throw (error.value ?? new Error(fallbackErrorMessage)) as Error
  }

  return data.value
}

/**
 * Requests one endpoint where only success/failure matters (response body may be empty).
 *
 * Accepted input examples:
 * - `endpoint = 'experiments/barcodes/'`, options = `{ method: 'POST', body: { experiment_id: 21, prefix: 'A001', number_of_plates: 4, sides: ['North'] } }`
 * - `endpoint = 'barcodespecifications/7/'`, options = `{ method: 'DELETE' }`
 */
export const requestApiVoid = async (
  endpoint: string,
  options: ApiRequestOptions,
  fallbackErrorMessage: string,
): Promise<void> => {
  const { error } = await useAPI<unknown>(endpoint, options)

  if (error.value) {
    throw (error.value ?? new Error(fallbackErrorMessage)) as Error
  }
}
