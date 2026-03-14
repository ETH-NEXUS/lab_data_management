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
    'text-white hover:text-green-900 border border-green-900 hover:border-green-600 bg-green-900 hover:bg-green-600',
  secondary:
    'text-slate-700 hover:text-slate-900 border border-slate-300 hover:border-slate-400 bg-white hover:bg-slate-100',
  info: 'text-white hover:text-blue-900 border border-blue-700 hover:border-blue-500 bg-blue-700 hover:bg-blue-500',
  danger: 'text-white hover:text-red-900 border border-red-700 hover:border-red-500 bg-red-700 hover:bg-red-500',
  accent:
    'text-slate-900 hover:text-slate-900 border border-lime-500 hover:border-lime-400 bg-[#BEF264] hover:bg-lime-300',
  outline:
    'uppercase tracking-wide text-slate-700 border border-transparent bg-transparent hover:bg-slate-100 hover:border-slate-200',
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
    'inline-flex items-center justify-center font-medium rounded-full transition duration-200 disabled:opacity-60 disabled:cursor-not-allowed',
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
