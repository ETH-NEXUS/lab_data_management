export default defineNuxtPlugin(() => {
  const authStore = useAuthStore()
  const authErrorStatus = useState<number | null>('authErrorStatus', () => null)

  watch(authErrorStatus, (status) => {
    if (status === 401 || status === 403) {
      authStore.clearUser()

      const route = useRoute()
      if (route.path !== '/login') {
        const next = encodeURIComponent(route.fullPath)
        navigateTo(`/login?next=${next}`)
      }
    }
  })
})