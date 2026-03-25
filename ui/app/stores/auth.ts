import { defineStore } from 'pinia'
import { ref, computed, watch } from 'vue'
import { useAPI } from '~~/composables/useAPI'

export type AuthUser = {
  id?: number
  username?: string
  email?: string
  first_name?: string
  last_name?: string
  groups?: string[]
  [key: string]: unknown
}

/**
 * -----------------------------
 * API ENDPOINTS
 * -----------------------------

 */
const LOGIN_URL = 'auth/login/'
const LOGOUT_URL = 'auth/logout/'
const ME_URL = 'auth/users/me/'

export const useAuthStore = defineStore('authStore', () => {
  const toast = useToast()

  const user = ref<AuthUser | null>(null)
  const isLoading = ref(false)
  const error = ref<string | null>(null)

  const isAuthenticated = computed(() => user.value != null)

  // global auth error status
  const authErrorStatus = useState<number | null>('authErrorStatus', () => null)

  // CSRF reset state
  const csrfToken = useState<string | null>('csrfToken', () => null)

  const clearError = () => {
    error.value = null
  }

  const setAuthStatusFromError = (err: any) => {
    const status = err?.response?.status
    if (status === 401 || status === 403) {
      authErrorStatus.value = status
    }
  }

  const showError = (title: string, description: string) => {
    toast.add({
      title,
      description,
      color: 'error',
      duration: 5000,
    })
  }

  /**
   * -----------------------------
   * Load current user
   * -----------------------------

   */
  const loadUser = async () => {
    isLoading.value = true
    error.value = null

    try {
      const { data, error: apiError } = await useAPI<AuthUser>(ME_URL, {
        method: 'GET',
      })

      if (apiError.value || !data.value) {
        user.value = null
        return null
      }

      user.value = data.value
      authErrorStatus.value = null
      return user.value
    } catch (err: any) {
      setAuthStatusFromError(err)
      user.value = null
      return null
    } finally {
      isLoading.value = false
    }
  }

  const clearUser = () => {
    user.value = null
  }

  /**
   * -----------------------------
   * Login
   * -----------------------------

   */
  const login = async (username: string, password: string) => {
    isLoading.value = true
    error.value = null

    try {
      const { error: apiError } = await useAPI(LOGIN_URL, {
        method: 'POST',
        body: { username, password },
      })

      if (apiError.value) {
        setAuthStatusFromError(apiError.value)
        user.value = null
        error.value = apiError.value.message ?? 'Login failed'
        showError('Login failed', error.value)
        return false
      }

      // reset CSRF cache
      csrfToken.value = null
      authErrorStatus.value = null

      await loadUser()
      return true
    } catch (err: any) {
      setAuthStatusFromError(err)
      user.value = null
      error.value = err instanceof Error ? err.message : 'Login failed'
      showError('Login failed', error.value)
      return false
    } finally {
      isLoading.value = false
    }
  }

  /**
   * -----------------------------
   * Logout
   * -----------------------------

   */
  const logout = async (redirectToLogin = true) => {
    isLoading.value = true
    error.value = null

    try {
      await useAPI(LOGOUT_URL, { method: 'GET' })
    } finally {
      user.value = null
      csrfToken.value = null
      authErrorStatus.value = null
      isLoading.value = false
    }

    if (redirectToLogin) {
      await navigateTo('/login')
    }
  }


  watch(authErrorStatus, (status: number | null) => {
    // If any request globally marks auth error (401/403), immediately clear user.
    if (status === 401 || status === 403) {
      user.value = null
    }
  })

  return {
    user,
    isAuthenticated,
    isLoading,
    error,
    clearError,
    loadUser,
    login,
    logout,
    clearUser,
  }
})
