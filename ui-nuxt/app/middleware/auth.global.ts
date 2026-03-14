export default defineNuxtRouteMiddleware(async (to) => {
  const authStore = useAuthStore()

  // routes that don't require login
  const normalizedPath = to.path.replace(/\/+$/, '') || '/'
  const isLoginRoute = normalizedPath === '/login' || /\/login$/.test(normalizedPath)
  if (isLoginRoute) return

  // if we already have a user, ok
  if (authStore.isAuthenticated) return

  // try to load session from backend (cookie-based)
  await authStore.loadUser()

  // still not authenticated → redirect to login with ?next=
  if (!authStore.isAuthenticated) {
    return navigateTo({
      path: '/login',
      query: { next: to.fullPath },
    })
  }
})
