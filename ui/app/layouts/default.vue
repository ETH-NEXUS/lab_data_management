<script setup lang="ts">
import AppDrawer from '~/components/layouts/AppDrawer.vue'
import AppHeader from '~/components/layouts/AppHeader.vue'

const route = useRoute()
const drawerWidth = 400

const leftDrawerOpen = ref(false)
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

    <div :class="shellClass">
      <AppDrawer :open="leftDrawerOpen" :width="drawerWidth" :top-offset="drawerTopOffset" @close="closeLeftDrawer" />

      <main :class="pageContainerClass">
        <NuxtPage :key="route.fullPath" />
      </main>
    </div>
  </div>
</template>
