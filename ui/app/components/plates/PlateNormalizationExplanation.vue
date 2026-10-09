<script setup lang="ts">
import { computed } from 'vue'
import PlateCalculationExplanation from '~/components/plates/PlateCalculationExplanation.vue'
import type { Plate } from '~/types/lab'
import type { NormalizedDataset } from '~/types/plateCalculations'
import { countWellsWithoutLog10, formatSummaryNumber, getWellTypeLog10Median } from '~/utils/plateDatasets'
import { computeZPrime } from '~/utils/plateStats'

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

/**
 * How many wells were left empty because their value is -1 or below, and how
 * many of them are controls: the medians of the controls are taken without
 * them, so the 0 and the 1 of the scale can be off (e.g. after a background
 * correction many P wells are below 0).
 */
const skippedWarning = (): string | undefined => {
  const { source, label, negativeType, positiveType } = props.dataset
  const count = countWellsWithoutLog10(props.plate, source, label)
  if (count === 0) return undefined

  const negativeCount = countWellsWithoutLog10(props.plate, source, label, negativeType)
  const positiveCount = countWellsWithoutLog10(props.plate, source, label, positiveType)
  const warning = t('plates.calculations.log10.skipped', { count })
  if (negativeCount === 0 && positiveCount === 0) return warning

  const controls = t('plates.calculations.normalization.skipped_controls', {
    negative: negativeType,
    positive: positiveType,
    negativeCount,
    positiveCount,
  })
  return `${warning} ${controls}`
}

/**
 * A warning if the controls overlap on the normalized values: their z' below 0
 * means the spread of their wells is larger than the distance between their
 * medians, e.g. when a wrong well type was chosen as a control. The z' does not
 * change with the scale, so it is the same for the %Inhibition and the %Activity.
 */
const overlapWarning = (): string | undefined => {
  const { label, negativeType, positiveType } = props.dataset
  const zPrime = computeZPrime(props.plate, label, positiveType, negativeType, props.timestampIndex)
  if (zPrime === null || zPrime >= 0) return undefined
  return t('plates.calculations.normalization.controls_overlap', {
    negative: negativeType,
    positive: positiveType,
    zPrime: formatSummaryNumber(zPrime),
  })
}

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
  }
  return {
    title: t(`plates.calculations.normalization.${kind}_title`, { label: props.dataset.label }),
    formulas: [
      t(`plates.calculations.normalization.${kind}_formula`, params),
      t(`plates.calculations.normalization.${kind}_formula_numbers`, params),
    ],
    reading: t(`plates.calculations.normalization.${kind}_reading`, params),
    // Both warnings, if there are any, e.g. "20 well(s) left empty ... The N1 and P controls are not clearly apart ..."
    warning: [skippedWarning(), overlapWarning()].filter(Boolean).join(' ') || undefined,
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
