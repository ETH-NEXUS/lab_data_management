<script setup lang="ts">
import { computed } from 'vue'
import PlateDatasetSummary from '~/components/plates/PlateDatasetSummary.vue'
import type { Plate } from '~/types/lab'

/**
 * What a heatmap shows, in plain words, below it: the title, which read the
 * numbers are of, what the values are, the formula (and the formula with the
 * numbers of the shown time point), the steps, how to read the values, and their
 * statistics. Raw data have no formula and no steps.
 *
 * Accepted props example:
 * - `{ title: 'log10: Lum1_log10', formulas: ['log10 value = log10(value of the well in Lum1)'],
 *     steps: ['Every well ...'], reading: '+1 means 10 times more signal ...', label: 'Lum1_log10' }`
 */
const props = withDefaults(
  defineProps<{
    plate: Plate
    label: string
    timestampIndex: number
    title: string
    // What the values are, e.g. for raw data
    description?: string
    formulas?: string[]
    steps?: string[]
    reading: string
    // e.g. how many wells were left empty
    warning?: string
  }>(),
  { description: undefined, formulas: () => [], steps: () => [], warning: undefined },
)

const { t } = useI18n()

// e.g. "Numbers of read 2 of 3 (2026-09-30T15:48:00Z)", or nothing without reads
const readCaption = computed(() => {
  const timestamps = props.plate.details.measurement_timestamps[props.label] ?? []
  const time = timestamps[props.timestampIndex]
  if (!time) return ''
  return t('plates.calculations.read_caption', {
    number: props.timestampIndex + 1,
    total: timestamps.length,
    time,
  })
})
</script>

<template>
  <div class="mt-3 space-y-1">
    <h4 class="font-medium text-slate-800">{{ props.title }}</h4>
    <p v-if="readCaption" class="text-xs text-slate-500">{{ readCaption }}</p>
    <p v-if="props.description" class="text-sm text-slate-700">{{ props.description }}</p>
    <template v-if="props.formulas.length > 0">
      <p class="text-sm text-slate-600">{{ t('plates.calculations.formula_caption') }}</p>
      <!-- Monospace with spaces, so "-" reads as a minus and not as a dash -->
      <div class="w-fit rounded-md bg-slate-100 px-3 py-2 font-mono text-sm text-slate-800 italic">
        <p v-for="formula in props.formulas" :key="formula">{{ formula }}</p>
      </div>
    </template>
    <template v-if="props.steps.length > 0">
      <p class="text-sm text-slate-600">{{ t('plates.calculations.steps_caption') }}</p>
      <ol class="list-decimal space-y-0.5 pl-6 text-sm text-slate-700">
        <li v-for="step in props.steps" :key="step">{{ step }}</li>
      </ol>
    </template>
    <p class="text-sm text-slate-700">
      <span class="font-medium">{{ t('plates.calculations.reading_caption') }}</span> {{ props.reading }}
    </p>
    <p v-if="props.warning" class="text-xs text-amber-700">{{ props.warning }}</p>
    <PlateDatasetSummary :plate="props.plate" :label="props.label" :timestamp-index="props.timestampIndex" />
  </div>
</template>
