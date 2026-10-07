<script setup lang="ts">
import { computed } from 'vue'
import type { ActivityDataset } from '~/types/backgroundCorrection'
import type { Plate } from '~/types/lab'
import { formatSummaryNumber, getWellTypeMedian } from '~/utils/plateDatasets'

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

const negativeMedian = computed(() =>
  getWellTypeMedian(props.plate, props.dataset.source, props.dataset.negativeType, props.timestampIndex),
)
const positiveMedian = computed(() =>
  getWellTypeMedian(props.plate, props.dataset.source, props.dataset.positiveType, props.timestampIndex),
)

const names = computed(() => ({
  negative: props.dataset.negativeType,
  positive: props.dataset.positiveType,
  source: props.dataset.source,
  negativeMedian: formatSummaryNumber(negativeMedian.value),
  positiveMedian: formatSummaryNumber(positiveMedian.value),
}))
</script>

<template>
  <div class="my-1 space-y-0.5 text-sm text-slate-700 italic">
    <p>{{ t('plates.calculations.activity_formula', names) }}</p>
    <p>{{ t('plates.calculations.activity_formula_numbers', names) }}</p>
    <p class="text-xs text-slate-500">{{ t('plates.calculations.activity_note', names) }}</p>
  </div>
</template>
