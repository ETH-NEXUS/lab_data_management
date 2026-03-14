<script setup lang="ts">
import AppDrawer from '~/components/layouts/AppDrawer.vue'
import AppHeader from '~/components/layouts/AppHeader.vue'

const route = useRoute()

const leftDrawerOpen = ref(false)
const toggleLeftDrawer = () => (leftDrawerOpen.value = !leftDrawerOpen.value)
const closeLeftDrawer = () => (leftDrawerOpen.value = false)
</script>

<template>
  <div class="min-h-screen bg-slate-50">
    <AppHeader :drawer-open="leftDrawerOpen" @toggle-drawer="toggleLeftDrawer" />

    <div class="flex min-h-[calc(100dvh-4rem)]">
      <AppDrawer :open="leftDrawerOpen" @close="closeLeftDrawer" />

      <main class="min-w-0 flex-1">
        <NuxtPage :key="route.fullPath" />
      </main>
    </div>
  </div>
</template>
