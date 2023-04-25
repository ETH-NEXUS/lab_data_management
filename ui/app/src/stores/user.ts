import {defineStore} from 'pinia'
import {api} from 'src/boot/axios'
import {ref} from 'vue'

interface Endpoints {
  loginCookie: string
  sessionLogin: string
  sessionLogout: string
  obtainToken: string
  refreshToken: string
  user: string
}

interface User {
  id: number
  email: string
  first_name: string
  last_name: string
  groups: Array<string>
}

interface ObtainTokenPayload {
  username: string
  password: string
}

interface UpdateTokenPayload {
  access: string
  refresh: string
}

export const useUserStore = defineStore('user', () => {
  const authenticated = ref<string | null>(localStorage.getItem('authenticated'))
  const jwt = ref<string | null>(localStorage.getItem('jwt'))
  const refreshJwt = ref<string | null>(localStorage.getItem('refreshJwt'))
  const user = ref<User | null>(JSON.parse(localStorage.getItem('user') || 'null'))
  const endpoints: Endpoints = {
    loginCookie: '/api/auth/cookie/',
    sessionLogin: '/api/auth/login/',
    sessionLogout: '/api/auth/logout/',
    obtainToken: '/api/auth/token/',
    refreshToken: '/api/auth/token/refresh/',
    user: '/api/auth/users/me/',
  }

  const _updateToken = (payload: UpdateTokenPayload) => {
    localStorage.setItem('authenticated', 'true')
    authenticated.value = 'true'
    if (payload.access) {
      localStorage.setItem('jwt', payload.access)
      jwt.value = payload.access
    }
    if (payload.refresh) {
      localStorage.setItem('refreshJwt', payload.refresh)
      refreshJwt.value = payload.refresh
    }
  }

  const _removeToken = () => {
    localStorage.removeItem('authenticated')
    localStorage.removeItem('jwt')
    localStorage.removeItem('refreshJwt')
    authenticated.value = null
    jwt.value = null
    refreshJwt.value = null
    user.value = null
  }

  const updateUserInfo = (payload: User) => {
    localStorage.setItem('user', JSON.stringify(payload))
    user.value = payload
  }

  const _removeUserInfo = () => {
    localStorage.removeItem('user')
    user.value = null
  }

  const obtainToken = async (payload: ObtainTokenPayload) => {
    try {
      const resp = await api.post(endpoints.obtainToken, payload)
      _updateToken(resp.data)
      await getUserInfo()
    } catch (err) {
      console.error(err)
    }
  }

  const sessionLogin = async (payload: ObtainTokenPayload) => {
    try {
      await api.get(endpoints.loginCookie)
      const resp = await api.post(endpoints.sessionLogin, payload)
      _updateToken(resp.data)
      await getUserInfo()
    } catch (err) {
      console.error(err)
    }
  }

  const sessionLogout = async () => {
    try {
      await api.get(endpoints.sessionLogout)
      removeToken()
    } catch (err) {
      console.error(err)
    }
  }

  const getUserInfo = async () => {
    try {
      const resp = await api.get(endpoints.user)
      updateUserInfo(resp.data)
    } catch (err) {
      console.error(err)
    }
  }

  const refreshToken = async () => {
    const payload = {
      refresh: refreshJwt.value,
    }
    try {
      const resp = await api.post(endpoints.refreshToken, payload)
      _updateToken(resp.data)
    } catch (err) {
      console.error(err)
    }
  }

  const removeToken = async () => {
    _removeToken()
    _removeUserInfo()
  }

  return {
    authenticated,
    jwt,
    refreshJwt,
    user,
    sessionLogin,
    sessionLogout,
    obtainToken,
    getUserInfo,
    refreshToken,
    removeToken,
  }
})
