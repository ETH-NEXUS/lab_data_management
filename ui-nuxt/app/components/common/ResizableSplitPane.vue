<script setup lang="ts">
import { computed, onBeforeUnmount, ref, watch } from 'vue'

type Props = {
  modelValue?: number
  minLeftPercent?: number
  maxLeftPercent?: number
  dividerWidth?: number
}

const props = withDefaults(defineProps<Props>(), {
  modelValue: 58,
  minLeftPercent: 25,
  maxLeftPercent: 75,
  dividerWidth: 4,
})

const emit = defineEmits<{
  (e: 'update:modelValue', value: number): void
}>()

const containerRef = ref<HTMLElement | null>(null)
const isDragging = ref(false)
const leftPercent = ref(props.modelValue)

/**
 * Clamps a split percentage so the divider always stays within safe bounds.
 *
 * Accepted data example:
 * - `15` with `minLeftPercent = 25` => clamped to `25`
 * - `82` with `maxLeftPercent = 75` => clamped to `75`
 *
 * Returned data example:
 * - `58`
 */
const clampPercent = (value: number): number => {
  const min = Math.min(props.minLeftPercent, props.maxLeftPercent)
  const max = Math.max(props.minLeftPercent, props.maxLeftPercent)
  return Math.min(max, Math.max(min, value))
}

watch(
  () => props.modelValue,
  (value) => {
    leftPercent.value = clampPercent(value)
  },
  { immediate: true },
)

/**
 * Applies a new split value locally and emits it to the parent page.
 *
 * Accepted data example:
 * - `63.2`
 *
 * Emitted data example:
 * - `update:modelValue` with `63.2`
 */
const applyPercent = (value: number): void => {
  const nextValue = clampPercent(value)
  leftPercent.value = nextValue
  emit('update:modelValue', nextValue)
}

/**
 * Converts an absolute pointer X coordinate into a left-panel percentage.
 *
 * Accepted data example:
 * - `clientX = 640`, container left = `120`, width = `1000` => `52`
 */
const updateFromClientX = (clientX: number): void => {
  const container = containerRef.value
  if (!container) return

  const rect = container.getBoundingClientRect()
  if (rect.width <= 0) return

  const rawPercent = ((clientX - rect.left) / rect.width) * 100
  applyPercent(rawPercent)
}

const onPointerMove = (event: PointerEvent): void => {
  if (!isDragging.value) return
  updateFromClientX(event.clientX)
}

const stopDragging = (): void => {
  if (!isDragging.value) return
  isDragging.value = false
  window.removeEventListener('pointermove', onPointerMove)
  window.removeEventListener('pointerup', stopDragging)
}

const startDragging = (event: PointerEvent): void => {
  event.preventDefault()
  isDragging.value = true
  window.addEventListener('pointermove', onPointerMove)
  window.addEventListener('pointerup', stopDragging)
}

/**
 * Keyboard accessibility for resizing:
 * - Left arrow reduces by 2%
 * - Right arrow increases by 2%
 * - Home snaps to min
 * - End snaps to max
 */
const onHandleKeyDown = (event: KeyboardEvent): void => {
  if (event.key === 'ArrowLeft') {
    event.preventDefault()
    applyPercent(leftPercent.value - 2)
    return
  }

  if (event.key === 'ArrowRight') {
    event.preventDefault()
    applyPercent(leftPercent.value + 2)
    return
  }

  if (event.key === 'Home') {
    event.preventDefault()
    applyPercent(props.minLeftPercent)
    return
  }

  if (event.key === 'End') {
    event.preventDefault()
    applyPercent(props.maxLeftPercent)
  }
}

onBeforeUnmount(() => {
  stopDragging()
})

const gridTemplateColumns = computed(() => `${leftPercent.value.toFixed(2)}% ${props.dividerWidth}px minmax(0, 1fr)`)
</script>

<template>
  <section ref="containerRef" class="grid h-dvh w-full" :style="{ gridTemplateColumns }">
    <div class="min-w-0">
      <slot name="left" />
    </div>

    <div
      role="separator"
      aria-orientation="vertical"
      tabindex="0"
      class="group relative cursor-col-resize bg-black/12 transition-colors hover:bg-black/25 focus:bg-black/25 focus:outline-none"
      @pointerdown="startDragging"
      @keydown="onHandleKeyDown"
    >
      <div class="pointer-events-none absolute inset-y-0 left-1/2 w-[1px] -translate-x-1/2 bg-black/35" />
      <div
        class="pointer-events-none absolute top-1/2 left-1/2 h-12 w-[2px] -translate-x-1/2 -translate-y-1/2 rounded-full bg-teal-900/30 transition-colors group-hover:bg-teal-900/60"
      />
    </div>

    <div class="min-w-0">
      <slot name="right" />
    </div>
  </section>
</template>
