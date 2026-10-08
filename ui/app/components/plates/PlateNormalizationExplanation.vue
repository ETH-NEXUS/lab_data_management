<script setup lang="ts">
import { computed } from 'vue'
import PlateCalculationExplanation from '~/components/plates/PlateCalculationExplanation.vue'
import type { Plate } from '~/types/lab'
import type { NormalizedDataset } from '~/types/plateCalculations'
import { formatSummaryNumber, getWellTypeLog10Median } from '~/utils/plateDatasets'

/**
 * How the %Inhibition or %Activity of the normalization was calculated, with the
 * medians of the controls of the shown time point, e.g. for `Lum1_activity_N_P`.
 */
const props = defineProps<{
  plate: Plate
  dataset: NormalizedDataset
  timestampIndex: number
}>()

const { t } = useI18n()

const texts = computed(() => {
  const { source, negativeType, positiveType, kind } = props.dataset
  const log10Median = (wellType: string) =>
    formatSummaryNumber(getWellTypeLog10Median(props.plate, source, wellType, props.timestampIndex))
  const params = {
    source,
    negative: negativeType,
    positive: positiveType,
    negativeMedian: log10Median(negativeType),
    positiveMedian: log10Median(positiveType),
    inhibitionLabel: props.dataset.inhibitionLabel,
    activityLabel: props.dataset.activityLabel,
  }
  // The same steps lead to both; %Activity adds the last one
  const steps = [
    t('plates.calculations.normalization.step_log', params),
    t('plates.calculations.normalization.step_negative', params),
    t('plates.calculations.normalization.step_positive', params),
    t('plates.calculations.normalization.step_inhibition', params),
  ]
  if (kind === 'activity') {
    steps.push(t('plates.calculations.normalization.step_activity', params))
  }
  return {
    title: t(`plates.calculations.normalization.${kind}_title`, { label: props.dataset.label }),
    formulas: [
      t(`plates.calculations.normalization.${kind}_formula`, params),
      t(`plates.calculations.normalization.${kind}_formula_numbers`, params),
    ],
    steps,
    reading: t(`plates.calculations.normalization.${kind}_reading`, params),
  }
})
</script>

<template>
  <PlateCalculationExplanation
    :plate="props.plate"
    :label="props.dataset.label"
    :timestamp-index="props.timestampIndex"
    v-bind="texts"
  />
</template>
