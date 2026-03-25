export const useCoreStore = defineStore('coreStore', () => {
  const colorMode = shallowRef<ReturnType<typeof useColorMode> | null>(null)

  const isDark = computed<boolean>(() => colorMode.value?.value === 'dark')

  const initColorMode = (): void => {
    if (import.meta.client && !colorMode.value) {
      colorMode.value = useColorMode()
    }
  }

  const toggleDarkMode = (value: boolean): void => {
    if (colorMode.value) {
      colorMode.value.preference = value ? 'dark' : 'light'
    }
  }

  const initFavicon = (): void => {
    if (!import.meta.client) return

    // set favicon to match system dark mode (browser UI), independent of app dark mode
    const isSystemDark = usePreferredDark()
    const favicon = computed<string>(() =>
      isSystemDark.value ? '/favicon/nexus_logo_dark_mode.png' : '/favicon/nexus_logo.png',
    )

    useFavicon(favicon, { rel: 'icon' })
  }

  return {
    colorMode,
    isDark,
    initColorMode,
    toggleDarkMode,
    initFavicon,
  }
})