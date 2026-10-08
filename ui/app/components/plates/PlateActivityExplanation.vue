<script setup lang="ts">
import { computed } from 'vue'
import PlateDatasetSummary from '~/components/plates/PlateDatasetSummary.vue'
import type { Plate } from '~/types/lab'
import type { ActivityDataset } from '~/types/plateCalculations'
import { formatSummaryNumber, getWellTypeStatistic } from '~/utils/plateDatasets'

/**
 * How the %Activity shown in the heatmap was calculated, with the numbers of
 * the shown time point put in, e.g. "= 100 · (value − 673) / (2958 − 673)".
 */
const props = defineProps<{
  plate: Plate
  dataset: ActivityDataset
  timestampIndex: number
}>()

const { t } = useI18n()

const medianOf = (wellType: string) =>
  getWellTypeStatistic(props.plate, props.dataset.source, wellType, 'median', props.timestampIndex)

// The parameters of the texts, e.g. { negative: 'N', negativeMedian: '2958', ... }
const formulaParams = computed(() => ({
  negative: props.dataset.negativeType,
  positive: props.dataset.positiveType,
  source: props.dataset.source,
  negativeMedian: formatSummaryNumber(medianOf(props.dataset.negativeType)),
  positiveMedian: formatSummaryNumber(medianOf(props.dataset.positiveType)),
}))
</script>

<template>
  <div class="mt-3">
    <h4 class="mb-1 font-medium text-slate-800">
      {{ t('plates.calculations.activity_title', { label: props.dataset.label }) }}
    </h4>
    <p class="text-sm text-slate-600">{{ t('plates.calculations.formula_caption') }}</p>
    <div class="my-1 space-y-0.5 text-sm text-slate-700 italic">
      <p>{{ t('plates.calculations.activity_formula', formulaParams) }}</p>
      <p>{{ t('plates.calculations.activity_formula_numbers', formulaParams) }}</p>
      <p class="text-xs text-slate-500">{{ t('plates.calculations.activity_note', formulaParams) }}</p>
    </div>
    <PlateDatasetSummary :plate="props.plate" :label="props.dataset.label" :timestamp-index="props.timestampIndex" />
  </div>
</template>
