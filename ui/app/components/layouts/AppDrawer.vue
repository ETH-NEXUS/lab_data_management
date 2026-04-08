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
  topOffset: 96,
})

const emit = defineEmits<{
  (e: 'close'): void
}>()

const onBackdrop = () => emit('close')
const onCloseClick = () => emit('close')

const route = useRoute()

/**
 * Determines if route navigation should close the drawer automatically.
 *
 * Returned value examples:
 * - `true`: mobile viewport (`< 768px`) where drawer should collapse after navigation.
 * - `false`: desktop viewport where the drawer state should stay as-is.
 */
const shouldAutoCloseOnRouteChange = (): boolean => {
  if (!import.meta.client) {
    return false
  }

  return !window.matchMedia('(min-width: 768px)').matches
}

watch(
  () => route.fullPath,
  () => {
    if (props.open && shouldAutoCloseOnRouteChange()) {
      emit('close')
    }
  },
)

/**
 * Drawer style variables used by responsive positioning classes.
 *
 * Accepted data example:
 * - `{ width: 400, topOffset: 96 }`
 *
 * Returned data example:
 * - `{ width: '400px', '--drawer-top-offset': '96px' }`
 */
const drawerStyle = computed<Record<string, string>>(() => ({
  width: `${props.width}px`,
  '--drawer-top-offset': `${props.topOffset}px`,
}))
</script>

<template>
  <!-- Backdrop (mobile) -->
  <div
    v-if="props.open"
    class="fixed inset-x-0 bottom-0 z-[70] bg-slate-900/30 md:hidden"
    :style="{ top: `${props.topOffset}px` }"
    @click="onBackdrop"
  />

  <!-- Drawer -->
  <aside
    id="app-navigation-drawer"
    class="fixed top-[var(--drawer-top-offset)] left-0 z-[75] h-[calc(100dvh-var(--drawer-top-offset))] shrink-0 border-r border-[var(--app-border)] bg-[var(--app-surface)]/95 shadow-[8px_0_30px_rgba(15,23,42,0.08)] backdrop-blur-xl transition-transform duration-200 ease-out md:static md:top-auto md:z-[40] md:h-auto md:bg-[var(--app-surface)] md:shadow-none md:backdrop-blur-none"
    :style="drawerStyle"
    :class="[props.open ? 'translate-x-0' : '-translate-x-full md:hidden']"
  >
    <div class="flex h-full flex-col md:h-auto">
      <div class="flex items-center justify-end border-b border-[var(--app-border)] p-2">
        <button
          type="button"
          aria-label="Close navigation"
          class="inline-flex cursor-pointer items-center justify-center rounded-md border border-blue-100 bg-white p-1.5 text-blue-700 transition hover:border-blue-200 hover:bg-blue-50 hover:text-blue-800"
          @click="onCloseClick"
        >
          <UIcon name="i-heroicons-x-mark" class="h-5 w-5" />
        </button>
      </div>
      <div class="min-h-0 flex-1 overflow-y-auto md:overflow-visible">
        <NavigationTree />
      </div>
    </div>
  </aside>
</template>
