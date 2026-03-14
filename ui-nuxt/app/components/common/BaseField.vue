<script setup lang="ts">
type Props = {
  modelValue: string
  label: string
  placeholder?: string
  type?: string
  autocomplete?: string
  fieldClass?: string // Wrapper around the whole field (label + input)
  labelClass?: string
  inputClass?: string
  hasRightSlot?: boolean //  Adds extra right padding for a right slot (eye icon)
}

const props = withDefaults(defineProps<Props>(), {
  placeholder: '',
  type: 'text',
  autocomplete: '',
  fieldClass: '',
  labelClass: 'block pl-4 mb-1 text-sm font-medium',
  inputClass:
    'w-full px-4 py-3 bg-white outline-none ring-offset-0 focus:ring-2 focus:ring-lime-500 shadow rounded-full',
  hasRightSlot: false
})

const emit = defineEmits<{
  (e: 'update:modelValue', value: string): void
}>()

const onInput = (e: Event) => {
  const target = e.target as HTMLInputElement
  emit('update:modelValue', target.value)
}
</script>

<template>
  <div :class="props.fieldClass">
    <label :class="props.labelClass">{{ props.label }}</label>

    <div class="relative">
      <input
        :value="props.modelValue"
        @input="onInput"
        :type="props.type"
        :autocomplete="props.autocomplete"
        :placeholder="props.placeholder"
        :class="[props.inputClass, props.hasRightSlot ? 'pr-12' : '']"
      />

      <!-- Right-side slot (eye icon, clear button, etc.) -->
      <div v-if="$slots.right" class="absolute top-1/2 right-0 -translate-y-1/2 mr-4">
        <slot name="right" />
      </div>
    </div>
  </div>
</template>