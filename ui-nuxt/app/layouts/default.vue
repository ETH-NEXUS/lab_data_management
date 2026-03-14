<script setup lang="ts">
import AppDrawer from '~/components/layouts/AppDrawer.vue'
import AppHeader from '~/components/layouts/AppHeader.vue'

const route = useRoute()
const drawerWidth = 400

const leftDrawerOpen = ref(false)
const toggleLeftDrawer = () => (leftDrawerOpen.value = !leftDrawerOpen.value)
const closeLeftDrawer = () => (leftDrawerOpen.value = false)

/**
 * Plate page uses a focused work layout (no top app header), closer to the legacy UI.
 *
 * Route examples that return `true`:
 * - `/plates/ABC001`
 * - `/en/plates/ABC001`
 */
const hideHeader = computed(() => /\/plates\//.test(route.path))
const drawerTopOffset = computed(() => (hideHeader.value ? 0 : 64))
const pageContainerClass = computed(() =>
  hideHeader.value ? 'min-h-dvh min-w-0 md:pl-[400px]' : 'min-h-[calc(100dvh-4rem)] min-w-0 md:pl-[400px]',
)
</script>

<template>
  <div class="relative min-h-screen overflow-hidden">
    <img
      src="/assets/double-waves.png"
      alt=""
      class="pointer-events-none fixed inset-0 -z-10 h-full w-full object-cover"
    />

    <AppHeader v-if="!hideHeader" :drawer-open="leftDrawerOpen" @toggle-drawer="toggleLeftDrawer" />

    <AppDrawer :open="leftDrawerOpen" :width="drawerWidth" :top-offset="drawerTopOffset" @close="closeLeftDrawer" />

    <main :class="pageContainerClass">
      <NuxtPage :key="route.fullPath" />
    </main>
  </div>
</template>
