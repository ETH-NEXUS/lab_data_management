<script setup lang="ts">
import PlateDatasetSummary from '~/components/plates/PlateDatasetSummary.vue'
import type { Plate } from '~/types/lab'

/**
 * How a calculated measurement was calculated, in plain words, below its heatmap:
 * the formula (and the formula with the numbers of the shown time point), the
 * steps, how to read the values, and their statistics.
 *
 * Accepted props example:
 * - `{ title: 'log10: Lum1_log10', formulas: ['log10 value = log10(value of the well in Lum1)'],
 *     steps: ['Every well ...'], reading: '+1 means 10 times more signal ...', label: 'Lum1_log10' }`
 */
const props = defineProps<{
  plate: Plate
  label: string
  timestampIndex: number
  title: string
  formulas: string[]
  steps: string[]
  reading: string
  // e.g. how many wells were left empty
  warning?: string
}>()

const { t } = useI18n()
</script>

<template>
  <div class="mt-3 space-y-1">
    <h4 class="font-medium text-slate-800">{{ props.title }}</h4>
    <p class="text-sm text-slate-600">{{ t('plates.calculations.formula_caption') }}</p>
    <!-- Monospace with spaces, so "-" reads as a minus and not as a dash -->
    <div class="w-fit rounded-md bg-slate-100 px-3 py-2 font-mono text-sm text-slate-800 italic">
      <p v-for="formula in props.formulas" :key="formula">{{ formula }}</p>
    </div>
    <p class="text-sm text-slate-600">{{ t('plates.calculations.steps_caption') }}</p>
    <ol class="list-decimal space-y-0.5 pl-6 text-sm text-slate-700">
      <li v-for="step in props.steps" :key="step">{{ step }}</li>
    </ol>
    <p class="text-sm text-slate-700">
      <span class="font-medium">{{ t('plates.calculations.reading_caption') }}</span> {{ props.reading }}
    </p>
    <p v-if="props.warning" class="text-xs text-amber-700">{{ props.warning }}</p>
    <PlateDatasetSummary :plate="props.plate" :label="props.label" :timestamp-index="props.timestampIndex" />
  </div>
</template>
