
import { useRuntimeConfig } from '#app'

type HttpMethod =
  | 'GET'
  | 'POST'
  | 'PUT'
  | 'DELETE'
  | 'PATCH'
  | 'HEAD'
  | 'OPTIONS'

export interface APIOptions {
  method?: HttpMethod
  headers?: Record<string, string>
  body?: any
  params?: Record<string, string>
  responseType?: 'json' | 'text' | 'blob' | 'arrayBuffer'
}

interface APIResponse<T> {
  data: { value: T | null }
  error: { value: Error | null }
}

/**
 * Reads CSRF token from cookie.
 */
function getCsrfTokenFromCookie(): string | null {
  if (typeof document === 'undefined') return null
  const token = document.cookie.match(/(?:^|;\s*)csrftoken=([^;]+)/)?.[1]
  return token ? decodeURIComponent(token) : null
}

/**
 * Ensures CSRF cookie exists before unsafe requests.
 */
let csrfInitPromise: Promise<void> | null = null

async function ensureCsrfCookie(apiBase: string) {
  if (typeof window === 'undefined') return
  if (getCsrfTokenFromCookie()) return

  if (!csrfInitPromise) {
    csrfInitPromise = $fetch(`${apiBase}/auth/cookie/`, {
      method: 'GET',
      credentials: 'include',
    })
      .then(() => {})
      .finally(() => {
        csrfInitPromise = null
      })
  }

  await csrfInitPromise
}

/**
 * Main API composable
 */
export const useAPI = async <T = any>(
  endpoint: string,
  options: APIOptions = {}
): Promise<APIResponse<T>> => {
  const config = useRuntimeConfig()
  const toast = useToast()

  const method = (options.method ?? 'GET').toUpperCase() as HttpMethod
  const isUnsafe = ['POST', 'PUT', 'DELETE', 'PATCH'].includes(method)

  const publicBase = (config.public?.baseURL ?? '').replace(/\/$/, '')
  const apiBase = publicBase ? `${publicBase}/api` : '/api'

  if (isUnsafe) {
    await ensureCsrfCookie(apiBase)
  }

  const headers: Record<string, string> = {
    ...(options.headers ?? {}),
  }

  const isFormData =
    typeof FormData !== 'undefined' && options.body instanceof FormData

  if (!isFormData && !headers['Content-Type']) {
    headers['Content-Type'] = 'application/json'
  }

  if (isUnsafe) {
    const token = getCsrfTokenFromCookie()
    if (token) {
      headers['X-CSRFToken'] = token
    }
  }

  const cleanEndpoint = endpoint.replace(/^\//, '')
  let url = `${apiBase}/${cleanEndpoint}`

  if (options.params) {
    const queryParams = new URLSearchParams()
    Object.entries(options.params).forEach(([key, value]) => {
      queryParams.append(key, value)
    })
    const queryString = queryParams.toString()
    if (queryString) {
      url += `?${queryString}`
    }
  }

  try {
    const response = await $fetch<T>(url, {
      method,
      headers,
      body: options.body,
      credentials: 'include',
      responseType: options.responseType,
    })

    return {
      data: { value: response },
      error: { value: null },
    }
  } catch (error: any) {
    const status = error?.response?.status as number | undefined
    const message = error?.message ?? 'Unknown error'

    // Don’t spam users for auth errors — middleware will redirect
    const isAuthError = status === 401 || status === 403

    if (!isAuthError) {
      toast.add({
        title: status ? `Error ${status}` : 'Network error',
        description: message,
        color: 'error',
        duration: 5000,
      })
    }

    return {
      data: { value: null },
      error: { value: error },
    }
  }
}
