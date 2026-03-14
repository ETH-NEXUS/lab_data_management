<script setup lang="ts">
import { computed, onBeforeUnmount, onMounted, ref } from 'vue'
import { useQueryClient } from '@tanstack/vue-query'
import { useAPI } from '~/composables/useAPI'
import { useAuthStore } from '~/stores/auth'
import { PROJECTS_QUERY_KEY } from '~/types/projects'

type Props = {
  transparent?: boolean
  drawerOpen?: boolean
}

const props = withDefaults(defineProps<Props>(), {
  transparent: false,
  drawerOpen: false,
})

const emit = defineEmits<{
  (e: 'toggle-drawer'): void
}>()

const authStore = useAuthStore()
const queryClient = useQueryClient()
const toast = useToast()
const { t } = useI18n()
const runtimeConfig = useRuntimeConfig()

const scrolled = ref(false)
const version = ref('N/A')
const isRefreshing = ref(false)

const onScroll = () => {
  scrolled.value = window.scrollY > 8
}

const loadVersion = async () => {
  const { data, error } = await useAPI<{ version?: string }>('version/', { method: 'GET' })
  if (!error.value) {
    version.value = data.value?.version ?? 'N/A'
  }
}

onMounted(async () => {
  onScroll()
  window.addEventListener('scroll', onScroll, { passive: true })
  await loadVersion()
})

onBeforeUnmount(() => {
  window.removeEventListener('scroll', onScroll)
})

const displayName = computed(() => {
  const u = authStore.user
  if (!u) return ''
  const full = `${u.first_name ?? ''} ${u.last_name ?? ''}`.trim()
  return full || u.username || u.email || ''
})

const toggleDrawer = () => emit('toggle-drawer')
const goHome = () => navigateTo('/')
const goMessages = () => navigateTo('/messages')

const refresh = async () => {
  if (isRefreshing.value) return

  isRefreshing.value = true
  try {
    await useAPI('refresh/', { method: 'GET' })
    await queryClient.invalidateQueries({ queryKey: PROJECTS_QUERY_KEY })
  } finally {
    isRefreshing.value = false
  }
}

/**
 * Normalizes configured backend base URL for non-Nuxt endpoints.
 *
 * Accepted input examples:
 * - `http://localhost:5090`
 * - `http://localhost:5090/`
 * - `http://localhost:5090/api`
 *
 * Returned value examples:
 * - `http://localhost:5090`
 */
const normalizeBackendBaseUrl = (value: string): string => {
  return value.replace(/\/$/, '').replace(/\/api$/, '')
}

const backendBaseUrl = computed(() => {
  const explicitBackendUrl = String(runtimeConfig.public?.backendURL ?? '').trim()
  if (explicitBackendUrl) {
    return normalizeBackendBaseUrl(explicitBackendUrl)
  }

  const apiBaseUrl = String(runtimeConfig.public?.baseURL ?? '').trim()
  if (apiBaseUrl) {
    return normalizeBackendBaseUrl(apiBaseUrl)
  }

  return ''
})

/**
 * Builds an absolute backend URL for resources not rendered by Nuxt (`/docs`, `/notebook`).
 *
 * Behavior:
 * - If `NUXT_PUBLIC_API_URL` is set, use that backend origin.
 * - Otherwise, use current origin where dev proxy paths (`/docs`, `/notebook`) are expected.
 *
 * Returned URL examples:
 * - with `backendURL = 'http://localhost:5090'` + `path = '/docs/'` => `http://localhost:5090/docs/`
 * - with empty `baseURL` + `path = '/notebook/'` => `http://localhost:8091/notebook/`
 */
const toBackendUrl = (path: string): string => {
  const normalizedPath = path.startsWith('/') ? path : `/${path}`
  if (backendBaseUrl.value) {
    return `${backendBaseUrl.value}${normalizedPath}`
  }

  return `${window.location.protocol}//${window.location.host}${normalizedPath}`
}

const openNotebook = () => {
  window.open(toBackendUrl('/notebook/'), '_blank', 'noreferrer')
}

const openDocsPage = () => {
  window.open(toBackendUrl('/docs/'), '_blank', 'noreferrer')
}

const logout = async () => {
  await authStore.logout(false)
  toast.add({
    title: t('app.header.logout'),
    description: 'Logged out',
    color: 'success',
    duration: 2000,
  })
  await navigateTo('/login')
}

const dropdownItems = computed(() => [
  [{ label: `Version ${version.value}`, icon: 'i-heroicons-map-pin' }],
  [
    { label: t('app.header.notebook'), icon: 'i-heroicons-pencil-square', onSelect: openNotebook },
    { label: t('app.header.help'), icon: 'i-heroicons-question-mark-circle', onSelect: openDocsPage },
  ],
  [{ label: t('app.header.logout'), icon: 'i-heroicons-arrow-right-on-rectangle', onSelect: logout }],
])
</script>

<template>
  <header
    class="sticky top-0 z-50 w-full border-b transition-all duration-300"
    :class="[
      props.transparent
        ? scrolled
          ? 'border-black/5 bg-white/70 backdrop-blur'
          : 'border-transparent bg-transparent'
        : 'border-black/5 bg-white',
    ]"
  >
    <div class="mx-auto h-16 w-full px-4">
      <div class="flex h-full items-center justify-between">
        <div class="flex items-center gap-2">
          <button
            class="inline-flex h-9 w-9 cursor-pointer items-center justify-center rounded-md transition hover:bg-teal-900/10"
            aria-label="Toggle navigation"
            aria-controls="app-navigation-drawer"
            :aria-expanded="props.drawerOpen ? 'true' : 'false'"
            type="button"
            @click="toggleDrawer"
          >
            <UIcon name="i-heroicons-bars-3" class="h-5 w-5" />
          </button>

          <button
            class="cursor-pointer text-lg font-semibold text-teal-900 transition hover:text-teal-700"
            type="button"
            @click="goHome"
          >
            {{ t('app.header.title') }}
          </button>

          <button
            class="inline-flex h-9 w-9 cursor-pointer items-center justify-center rounded-md transition hover:bg-teal-900/10 disabled:cursor-not-allowed"
            :aria-label="t('app.header.refresh')"
            :disabled="isRefreshing"
            type="button"
            @click="refresh"
          >
            <UIcon name="i-heroicons-arrow-path" :class="['h-5 w-5', isRefreshing && 'animate-spin']" />
          </button>
        </div>

        <div class="flex items-center gap-3">
          <span v-if="authStore.isAuthenticated" class="hidden text-sm font-medium text-teal-900 md:inline">
            {{ displayName }}
          </span>

          <button
            class="inline-flex h-9 w-9 cursor-pointer items-center justify-center rounded-md transition hover:bg-teal-900/10"
            :aria-label="t('app.header.messages')"
            type="button"
            @click="goMessages"
          >
            <UIcon name="i-heroicons-flag" class="h-5 w-5" />
          </button>

          <UDropdownMenu
            v-if="authStore.isAuthenticated"
            :items="dropdownItems"
            :content="{ side: 'bottom', align: 'end' }"
            :ui="{ content: 'z-[9999]', item: 'cursor-pointer', label: 'cursor-pointer' }"
          >
            <button
              class="inline-flex h-9 w-9 cursor-pointer items-center justify-center rounded-md transition hover:bg-teal-900/10"
              aria-label="User menu"
              type="button"
            >
              <UIcon name="i-heroicons-user" class="h-5 w-5" />
            </button>
          </UDropdownMenu>
        </div>
      </div>
    </div>
  </header>
</template>
