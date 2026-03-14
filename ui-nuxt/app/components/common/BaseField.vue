<script setup lang="ts">
type Props = {
  modelValue: string
  label: string
  placeholder?: string
  type?: string
  autocomplete?: string
  autofocus?: boolean
  multiline?: boolean
  rows?: number
  fieldClass?: string // Wrapper around the whole field (label + input)
  labelClass?: string
  inputClass?: string
  hasRightSlot?: boolean //  Adds extra right padding for a right slot (eye icon)
}

const props = withDefaults(defineProps<Props>(), {
  placeholder: '',
  type: 'text',
  autocomplete: '',
  autofocus: false,
  multiline: false,
  rows: 5,
  fieldClass: '',
  labelClass: 'block pl-4 mb-1 text-sm font-medium',
  inputClass:
    'w-full px-4 py-3 bg-white outline-none ring-offset-0 focus:ring-2 focus:ring-lime-500 shadow rounded-full',
  hasRightSlot: false,
})

const emit = defineEmits<{
  (e: 'update:modelValue', value: string): void
}>()

const onInput = (e: Event) => {
  const target = e.target as HTMLInputElement | HTMLTextAreaElement
  emit('update:modelValue', target.value)
}
</script>

<template>
  <div :class="props.fieldClass">
    <label :class="props.labelClass">{{ props.label }}</label>

    <div class="relative">
      <textarea
        v-if="props.multiline"
        :value="props.modelValue"
        :rows="props.rows"
        :placeholder="props.placeholder"
        :class="[props.inputClass, props.hasRightSlot ? 'pr-12' : '']"
        :autofocus="props.autofocus"
        @input="onInput"
      />

      <input
        v-else
        :value="props.modelValue"
        :type="props.type"
        :autocomplete="props.autocomplete"
        :placeholder="props.placeholder"
        :autofocus="props.autofocus"
        :class="[props.inputClass, props.hasRightSlot ? 'pr-12' : '']"
        @input="onInput"
      />

      <!-- Right-side slot (eye icon, clear button, etc.) -->
      <div v-if="$slots.right" class="absolute top-1/2 right-0 mr-4 -translate-y-1/2">
        <slot name="right" />
      </div>
    </div>
  </div>
</template>
