<script setup lang="ts">
import AppDrawer from '~/components/layouts/AppDrawer.vue'
import AppHeader from '~/components/layouts/AppHeader.vue'

const route = useRoute()
const drawerWidth = 400

const leftDrawerOpen = ref(false)
const toggleLeftDrawer = () => (leftDrawerOpen.value = !leftDrawerOpen.value)
const closeLeftDrawer = () => (leftDrawerOpen.value = false)
const drawerTopOffset = 0
const pageContainerClass = 'min-h-[calc(100dvh-6rem)] min-w-0 md:pl-[400px]'
</script>

<template>
  <div class="relative min-h-screen overflow-hidden">
    <!-- Background image removed. Use lightweight gradients for depth. -->
    <div class="pointer-events-none fixed inset-0 -z-10 bg-[var(--app-bg)]" />
    <div
      class="pointer-events-none fixed -top-24 -left-16 -z-10 h-96 w-96 rounded-full bg-[var(--app-accent-soft)]/70 blur-3xl"
    />
    <div class="pointer-events-none fixed right-0 bottom-0 -z-10 h-96 w-96 rounded-full bg-white/70 blur-3xl" />

    <AppHeader :drawer-open="leftDrawerOpen" @toggle-drawer="toggleLeftDrawer" />

    <AppDrawer :open="leftDrawerOpen" :width="drawerWidth" :top-offset="drawerTopOffset" @close="closeLeftDrawer" />

    <main :class="pageContainerClass">
      <NuxtPage :key="route.fullPath" />
    </main>
  </div>
</template>
