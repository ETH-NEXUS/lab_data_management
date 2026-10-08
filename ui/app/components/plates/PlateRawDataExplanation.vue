<script setup lang="ts">
import { computed } from 'vue'
import PlateCalculationExplanation from '~/components/plates/PlateCalculationExplanation.vue'
import type { Plate } from '~/types/lab'
import type { PlateMeasurementSources } from '~/types/plateCalculations'

/**
 * What the heatmap of a measurement that was not calculated on this page shows,
 * e.g. `Lum1`: where its values come from (the file it was imported from, or the
 * measurement calculator of the experiment), what the colours mean and the
 * medians of the well types, e.g. of the R wells (the background).
 *
 * Accepted props example:
 * - `{ label: 'Lum1', sources: { Lum1: '093026-154654_RKS_300926_1.asc' }, ... }`
 */
const props = defineProps<{
  plate: Plate
  label: string
  timestampIndex: number
  // Not known yet while they are loaded: then the general text is shown
  sources: PlateMeasurementSources
}>()

const { t } = useI18n()

const texts = computed(() => {
  const label = props.label
  if (!(label in props.sources)) {
    return {
      title: t('plates.calculations.raw.title', { label }),
      description: t('plates.calculations.raw.description'),
    }
  }
  const file = props.sources[label]
  if (file) {
    return {
      title: t('plates.calculations.raw.title_imported', { label }),
      description: t('plates.calculations.raw.description_imported', { file }),
    }
  }
  return {
    title: t('plates.calculations.raw.title_calculator', { label }),
    description: t('plates.calculations.raw.description_calculator'),
  }
})
</script>

<template>
  <PlateCalculationExplanation
    :plate="props.plate"
    :label="props.label"
    :timestamp-index="props.timestampIndex"
    :title="texts.title"
    :description="texts.description"
    :reading="t('plates.calculations.raw.reading')"
  />
</template>
