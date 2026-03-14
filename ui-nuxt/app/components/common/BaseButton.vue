<script setup lang="ts">
type Props = {
  label: string
  onClick: () => void | Promise<void>
  loading?: boolean
  disabled?: boolean
  className?: string
}

const props = withDefaults(defineProps<Props>(), {
  loading: false,
  disabled: false,
  className:
    'inline-flex w-full py-3 px-6 items-center justify-center text-lg font-medium ' +
    'text-white hover:text-green-900 ' +
    'border border-green-900 hover:border-green-600 ' +
    'bg-green-900 hover:bg-green-600 ' +
    'rounded-full transition duration-200 ' +
    'disabled:opacity-60 disabled:cursor-not-allowed'
})

const isRunning = ref(false)

const handleClick = async () => {
  if (props.disabled || props.loading || isRunning.value) return

  try {
    isRunning.value = true
    await props.onClick()
  } finally {
    isRunning.value = false
  }
}

const showSpinner = computed(() => props.loading || isRunning.value)
</script>

<template>
  <button
    type="button"
    :class="props.className"
    :disabled="props.disabled || props.loading || isRunning"
    @click="handleClick"
  >
    <span v-if="showSpinner" class="mr-3">
      <span class="inline-block h-5 w-5 rounded-full border-2 border-white/70 border-t-transparent animate-spin" />
    </span>

    {{ props.label }}
  </button>
</template>