<script setup lang="ts">
/**
 * One choice out of a few, shown as radio buttons, e.g. the file name format of the M1000.
 *
 * Accepted props example:
 * - `{ label: 'File name format', choices: ['date_time_barcode', 'barcode_date_time'],
 *     choiceLabels: { date_time_barcode: '20240610-121212_demo_1.asc' }, modelValue: 'date_time_barcode' }`
 */
const props = defineProps<{
  label: string
  name: string
  choices: string[]
  choiceLabels?: Record<string, string>
  modelValue: string
}>()

const emit = defineEmits<{
  (e: 'update:modelValue', value: string): void
}>()
</script>

<template>
  <fieldset class="space-y-2">
    <legend class="mb-2 block pl-1 text-sm font-semibold tracking-[0.04em] text-slate-600">
      {{ props.label }}
    </legend>
    <label
      v-for="choice in props.choices"
      :key="`${props.name}-${choice}`"
      class="flex cursor-pointer items-center gap-2 pl-1 text-sm text-slate-700"
    >
      <input
        type="radio"
        :name="props.name"
        :value="choice"
        :checked="props.modelValue === choice"
        class="h-4 w-4 cursor-pointer text-blue-600 focus:ring-blue-300"
        @change="emit('update:modelValue', choice)"
      />
      <span>{{ props.choiceLabels?.[choice] ?? choice }}</span>
    </label>
  </fieldset>
</template>
