<script setup lang="ts">
import { computed } from 'vue'
import PlateCalculationExplanation from '~/components/plates/PlateCalculationExplanation.vue'
import type { Plate } from '~/types/lab'
import type { BackgroundCorrection } from '~/types/plateCalculations'
import {
  formatSummaryNumber,
  getCorrectionDataset,
  getLog10SourceLabel,
  getWellTypeStatistic,
} from '~/utils/plateDatasets'

/**
 * How a background correction was calculated, with the subtracted value of the
 * shown time point, e.g. for the correction `Lum1_bc_R_median` of `Lum1`.
 */
const props = defineProps<{
  plate: Plate
  correction: BackgroundCorrection
  timestampIndex: number
}>()

const { t } = useI18n()

const texts = computed(() => {
  const { source, referenceType, method } = props.correction
  // The value that was subtracted, e.g. the median of the R wells
  const background = getWellTypeStatistic(props.plate, source, referenceType, method, props.timestampIndex)
  const params = {
    source,
    reference: referenceType,
    method: t(`plates.background_correction.methods.${method}`),
    value: formatSummaryNumber(background),
  }
  // In log10 values, subtracting means dividing. A correction of a correction
  // (e.g. Lum1_log10_bc_R_median_bc_N_median) is of log10 values if the first one was.
  let firstSource = source
  let earlierCorrection = getCorrectionDataset(props.plate, firstSource)
  while (earlierCorrection) {
    firstSource = earlierCorrection.source
    earlierCorrection = getCorrectionDataset(props.plate, firstSource)
  }
  const isOfLog10 = getLog10SourceLabel(props.plate, firstSource) !== null
  return {
    title: t('plates.background_correction.heatmap_title', { label: props.correction.label }),
    formulas: [
      t('plates.background_correction.formula', params),
      t('plates.background_correction.formula_numbers', params),
    ],
    steps: [
      t('plates.background_correction.step_background', params),
      t('plates.background_correction.step_subtract', params),
      t('plates.background_correction.step_empty', params),
    ],
    reading: t(
      isOfLog10 ? 'plates.background_correction.reading_log10' : 'plates.background_correction.reading',
      params,
    ),
  }
})
</script>

<template>
  <PlateCalculationExplanation
    :plate="props.plate"
    :label="props.correction.label"
    :timestamp-index="props.timestampIndex"
    v-bind="texts"
  />
</template>
