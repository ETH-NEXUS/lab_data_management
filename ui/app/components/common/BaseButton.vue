<script setup lang="ts">
type ButtonVariant = 'primary' | 'secondary' | 'info' | 'danger' | 'accent' | 'outline'
type ButtonSize = 'lg' | 'sm'
type ButtonWidth = 'full' | 'auto'

type Props = {
  label: string
  onClick: () => void | Promise<void>
  loading?: boolean
  disabled?: boolean
  variant?: ButtonVariant
  size?: ButtonSize
  width?: ButtonWidth
  className?: string
}

const props = withDefaults(defineProps<Props>(), {
  loading: false,
  disabled: false,
  variant: 'primary',
  size: 'lg',
  width: 'full',
  className: '',
})

const isRunning = ref(false)

/**
 * Shared style presets for reusable button variants.
 *
 * Data examples:
 * - `{ variant: 'primary', size: 'lg', width: 'full' }`
 * - `{ variant: 'secondary', size: 'sm', width: 'auto' }`
 * - `{ variant: 'danger', size: 'sm', width: 'auto' }`
 * - `{ variant: 'accent', size: 'sm', width: 'auto' }`
 * - `{ variant: 'outline', size: 'sm', width: 'auto' }`
 */
const variantClassMap: Record<ButtonVariant, string> = {
  primary:
    'border border-[var(--app-accent)] bg-[var(--app-accent)] text-white hover:border-[var(--app-accent-hover)] hover:bg-[var(--app-accent-hover)]',
  secondary:
    'border border-[var(--app-border)] bg-[var(--app-surface)] text-[var(--app-text-primary)] hover:bg-[var(--app-surface-muted)]',
  info: 'border border-blue-700 bg-blue-700 text-white hover:border-blue-600 hover:bg-blue-600',
  danger: 'border border-red-700 bg-red-700 text-white hover:border-red-600 hover:bg-red-600',
  accent:
    'border border-amber-400 bg-amber-300 text-slate-900 hover:border-amber-300 hover:bg-amber-200',
  outline:
    'border border-transparent bg-transparent uppercase tracking-wide text-[var(--app-text-muted)] hover:border-[var(--app-border)] hover:bg-[var(--app-surface)]',
}

const sizeClassMap: Record<ButtonSize, string> = {
  lg: 'py-3 px-6 text-lg',
  sm: 'py-2 px-4 text-sm',
}

const widthClassMap: Record<ButtonWidth, string> = {
  full: 'w-full',
  auto: 'w-auto',
}

const handleClick = async () => {
  if (props.disabled || props.loading || isRunning.value) return

  try {
    isRunning.value = true
    await props.onClick()
  } finally {
    isRunning.value = false
  }
}

const buttonClass = computed(() => {
  return [
    'inline-flex items-center justify-center rounded-xl font-medium transition duration-200 disabled:cursor-not-allowed disabled:opacity-60',
    variantClassMap[props.variant],
    sizeClassMap[props.size],
    widthClassMap[props.width],
    props.className,
  ].join(' ')
})

const showSpinner = computed(() => props.loading || isRunning.value)
const spinnerClass = computed(() => (props.size === 'sm' ? 'h-4 w-4' : 'h-5 w-5'))
const spinnerColorClass = computed(() => {
  if (props.variant === 'secondary' || props.variant === 'outline') return 'border-slate-500/70'
  if (props.variant === 'accent') return 'border-slate-700/70'
  return 'border-white/70'
})
</script>

<template>
  <button
    type="button"
    :class="buttonClass"
    :disabled="props.disabled || props.loading || isRunning"
    @click="handleClick"
  >
    <span v-if="showSpinner" class="mr-3">
      <span
        :class="[
          'inline-block animate-spin rounded-full border-2 border-t-transparent',
          spinnerClass,
          spinnerColorClass,
        ]"
      />
    </span>

    {{ props.label }}
  </button>
</template>
