<script setup lang="ts">
import AppDrawer from '~/components/layouts/AppDrawer.vue'
import AppHeader from '~/components/layouts/AppHeader.vue'

const route = useRoute()
const drawerWidth = 400

/**
 * Tracks whether the navigation drawer is visible.
 *
 * State examples:
 * - `true`: navigation drawer is rendered on screen.
 * - `false`: navigation drawer is hidden and can be reopened from the page button.
 */
const leftDrawerOpen = ref(true)
const openLeftDrawer = () => (leftDrawerOpen.value = true)
const closeLeftDrawer = () => (leftDrawerOpen.value = false)
const drawerTopOffset = 96
const shellClass = 'min-h-[calc(100dvh-6rem)] md:flex md:items-stretch'
const pageContainerClass = 'min-w-0 flex-1 md:border-l md:border-[var(--app-border)]'
</script>

<template>
  <div class="relative min-h-screen">
    <!-- Background image removed. Use lightweight gradients for depth. -->
    <div class="pointer-events-none fixed inset-0 -z-10 bg-[var(--app-bg)]" />
    <div
      class="pointer-events-none fixed -top-24 -left-16 -z-10 h-96 w-96 rounded-full bg-[var(--app-accent-soft)]/70 blur-3xl"
    />
    <div class="pointer-events-none fixed right-0 bottom-0 -z-10 h-96 w-96 rounded-full bg-white/70 blur-3xl" />

    <AppHeader />

    <button
      v-if="!leftDrawerOpen"
      type="button"
      aria-controls="app-navigation-drawer"
      aria-label="Open navigation"
      class="fixed left-3 z-[80] inline-flex cursor-pointer items-center gap-2 rounded-full border border-blue-100 bg-white px-3 py-2 text-xs font-semibold tracking-[0.1em] text-blue-700 uppercase shadow-sm transition hover:border-blue-200 hover:bg-blue-50 hover:text-blue-800 sm:left-4"
      :style="{ top: `${drawerTopOffset + 12}px` }"
      @click="openLeftDrawer"
    >
      <UIcon name="i-heroicons-bars-3" class="h-4 w-4" />
      <span>Navigation</span>
    </button>

    <div :class="shellClass">
      <AppDrawer :open="leftDrawerOpen" :width="drawerWidth" :top-offset="drawerTopOffset" @close="closeLeftDrawer" />

      <main :class="pageContainerClass">
        <NuxtPage :key="route.fullPath" />
      </main>
    </div>
  </div>
</template>
