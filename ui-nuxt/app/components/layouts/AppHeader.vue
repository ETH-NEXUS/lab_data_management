<script setup lang="ts">
import { computed, onMounted, ref } from 'vue'
import { useQueryClient } from '@tanstack/vue-query'
import { useAPI } from '~/composables/useAPI'
import { useAuthStore } from '~/stores/auth'
import { PROJECTS_QUERY_KEY } from '~/types/projects'

type Props = {
  drawerOpen?: boolean
}

const props = withDefaults(defineProps<Props>(), {
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
const route = useRoute()

const version = ref('N/A')
const isRefreshing = ref(false)

const iconButtonClass =
  'inline-flex h-10 w-10 cursor-pointer items-center justify-center rounded-xl border border-white/30 text-white transition hover:border-white/70 hover:bg-white/20'
const iconButtonDisabledClass = 'disabled:cursor-not-allowed disabled:opacity-60'
const menuButtonClass =
  'inline-flex h-10 w-10 cursor-pointer items-center justify-center rounded-xl border border-white/40 text-white transition hover:border-white/80 hover:bg-white/20'

const navLinks = computed(() => [
  { label: 'Projects', to: '/projects' },
  { label: 'Experiments', to: '/experiments' },
  { label: t('app.header.messages'), to: '/messages' },
  { label: 'Management', to: '/management' },
])

/**
 * Checks if a navigation item is active for prefixed and nested routes.
 *
 * Accepted data examples:
 * - `{ routePath: '/projects', to: '/projects' }`
 * - `{ routePath: '/en/projects/42', to: '/projects' }`
 *
 * Returned examples:
 * - `true`
 * - `false`
 */
const isNavActive = (to: string): boolean => {
  return route.path === to || route.path.endsWith(to) || route.path.includes(`${to}/`)
}

const navLinkClass = (to: string): string => {
  if (isNavActive(to)) {
    return 'inline-flex rounded-full bg-white px-4 py-2 text-xs font-semibold tracking-[0.12em] text-blue-700 uppercase transition'
  }

  return 'inline-flex rounded-full border border-white/35 bg-white/12 px-4 py-2 text-xs font-semibold tracking-[0.12em] text-white uppercase transition hover:border-blue-100 hover:bg-white hover:text-blue-700'
}

const loadVersion = async () => {
  const { data, error } = await useAPI<{ version?: string }>('version/', { method: 'GET' })
  if (!error.value) {
    version.value = data.value?.version ?? 'N/A'
  }
}

onMounted(async () => {
  await loadVersion()
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
  <header class="sticky top-0 z-50 w-full shadow-md">
    <div class="border-b border-blue-400/70 bg-blue-600">
      <div class="mx-auto flex h-8 w-full items-center justify-center px-3">
        <UIcon name="i-heroicons-sparkles" class="mr-2 h-3.5 w-3.5 text-blue-50" />
        <p class="text-[11px] font-semibold tracking-[0.12em] text-blue-50 uppercase">Refresh to sync latest data.</p>
      </div>
    </div>

    <nav class="app-nav-surface border-b border-blue-300/60">
      <div class="mx-auto flex h-16 w-full items-center justify-between gap-3 px-3 sm:px-4 lg:px-6">
        <div class="flex items-center gap-2">
          <button
            :class="menuButtonClass"
            aria-label="Toggle navigation"
            aria-controls="app-navigation-drawer"
            :aria-expanded="props.drawerOpen ? 'true' : 'false'"
            type="button"
            @click="toggleDrawer"
          >
            <UIcon name="i-heroicons-bars-3" class="h-5 w-5" />
          </button>

          <button
            class="cursor-pointer text-sm font-semibold tracking-[0.08em] text-white uppercase transition hover:text-blue-50 sm:text-base"
            type="button"
            @click="goHome"
          >
            {{ t('app.header.title') }}
          </button>
        </div>

        <ul class="hidden items-center gap-2 xl:flex">
          <li v-for="link in navLinks" :key="link.to">
            <NuxtLink :to="link.to" :class="navLinkClass(link.to)">
              {{ link.label }}
            </NuxtLink>
          </li>
        </ul>

        <div class="flex items-center gap-2">
          <span
            class="hidden rounded-full bg-white/18 px-3 py-1 text-[11px] font-semibold tracking-[0.12em] text-blue-50 uppercase 2xl:inline-flex"
          >
            v{{ version }}
          </span>

          <button
            :class="[iconButtonClass, iconButtonDisabledClass]"
            :aria-label="t('app.header.refresh')"
            :disabled="isRefreshing"
            type="button"
            @click="refresh"
          >
            <UIcon name="i-heroicons-arrow-path" :class="['h-5 w-5', isRefreshing && 'animate-spin']" />
          </button>

          <button :class="iconButtonClass" :aria-label="t('app.header.messages')" type="button" @click="goMessages">
            <UIcon name="i-heroicons-flag" class="h-5 w-5" />
          </button>

          <UDropdownMenu
            v-if="authStore.isAuthenticated"
            :items="dropdownItems"
            :content="{ side: 'bottom', align: 'end' }"
            :ui="{ content: 'z-[9999]', item: 'cursor-pointer', label: 'cursor-pointer' }"
          >
            <button :class="iconButtonClass" aria-label="User menu" type="button">
              <UIcon name="i-heroicons-user" class="h-5 w-5" />
            </button>
          </UDropdownMenu>

          <span v-if="authStore.isAuthenticated" class="hidden text-sm font-medium text-blue-50 2xl:inline">
            {{ displayName }}
          </span>
        </div>
      </div>
    </nav>
  </header>
</template>
