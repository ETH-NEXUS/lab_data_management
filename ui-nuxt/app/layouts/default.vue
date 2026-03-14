<script setup lang="ts">
import AppDrawer from '~/components/layouts/AppDrawer.vue'
import AppHeader from '~/components/layouts/AppHeader.vue'

const route = useRoute()
const drawerWidth = 400

const leftDrawerOpen = ref(false)
const toggleLeftDrawer = () => (leftDrawerOpen.value = !leftDrawerOpen.value)
const closeLeftDrawer = () => (leftDrawerOpen.value = false)
</script>

<template>
  <div class="relative min-h-screen overflow-hidden">
    <img
      src="/assets/double-waves.png"
      alt=""
      class="pointer-events-none fixed inset-0 -z-10 h-full w-full object-cover"
    />

    <AppHeader :drawer-open="leftDrawerOpen" @toggle-drawer="toggleLeftDrawer" />

    <AppDrawer :open="leftDrawerOpen" :width="drawerWidth" @close="closeLeftDrawer" />

    <main class="min-h-[calc(100dvh-4rem)] min-w-0 md:pl-[400px]">
      <NuxtPage :key="route.fullPath" />
    </main>
  </div>
</template>
