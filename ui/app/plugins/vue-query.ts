import { VueQueryPlugin, QueryClient, type VueQueryPluginOptions } from '@tanstack/vue-query'

export default defineNuxtPlugin((nuxtApp) => {
  const queryClient = new QueryClient({
    defaultOptions: {
      queries: {
        staleTime: 1000 * 60 * 5, // 5 minutes
        gcTime: 1000 * 60 * 30, // 30 minutes (formerly cacheTime)
        retry: 1,
        refetchOnWindowFocus: false,
      },
    },
  })

  const options: VueQueryPluginOptions = {
    queryClient,
  }

  nuxtApp.vueApp.use(VueQueryPlugin, options)
})
