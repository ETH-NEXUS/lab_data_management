<script setup lang="ts">
import PlateCalculationExplanation from '~/components/plates/PlateCalculationExplanation.vue'
import type { Plate } from '~/types/lab'

/**
 * What the heatmap of a measurement that was not calculated on this page shows
 * (the raw data, e.g. `Lum1`): no formula, only what the values and colours mean
 * and the medians of the well types, e.g. of the R wells (the background).
 */
const props = defineProps<{
  plate: Plate
  label: string
  timestampIndex: number
}>()

const { t } = useI18n()
</script>

<template>
  <PlateCalculationExplanation
    :plate="props.plate"
    :label="props.label"
    :timestamp-index="props.timestampIndex"
    :title="t('plates.calculations.raw.title', { label: props.label })"
    :description="t('plates.calculations.raw.description')"
    :reading="t('plates.calculations.raw.reading')"
  />
</template>
