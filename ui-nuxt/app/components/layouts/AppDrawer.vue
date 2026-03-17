<script setup lang="ts">
import { computed, watch } from 'vue'
import NavigationTree from '~/components/navigation/NavigationTree.vue'

type Props = {
  open: boolean
  width?: number
  topOffset?: number
}

const props = withDefaults(defineProps<Props>(), {
  width: 400,
  topOffset: 64,
})

const emit = defineEmits<{
  (e: 'close'): void
}>()

const onBackdrop = () => emit('close')

const route = useRoute()

watch(
  () => route.fullPath,
  () => {
    if (props.open) {
      emit('close')
    }
  },
)

const drawerStyle = computed(() => ({
  width: `${props.width}px`,
  top: `${props.topOffset}px`,
  height: `calc(100dvh - ${props.topOffset}px)`,
}))
</script>

<template>
  <!-- Backdrop (mobile) -->
  <div
    v-if="props.open"
    class="fixed inset-x-0 bottom-0 z-30 bg-slate-900/30 md:hidden"
    :style="{ top: `${props.topOffset}px` }"
    @click="onBackdrop"
  />

  <!-- Drawer -->
  <aside
    id="app-navigation-drawer"
    class="fixed left-0 z-40 shrink-0 border-r border-[var(--app-border)] bg-[var(--app-surface)]/95 shadow-[8px_0_30px_rgba(15,23,42,0.08)] backdrop-blur-xl transition-transform duration-200 ease-out md:translate-x-0"
    :style="drawerStyle"
    :class="[props.open ? 'translate-x-0' : '-translate-x-full']"
  >
    <div class="h-full overflow-auto">
      <NavigationTree />
    </div>
  </aside>
</template>
