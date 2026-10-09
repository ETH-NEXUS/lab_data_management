<script setup lang="ts">
import { computed } from 'vue'
import PlateCalculationExplanation from '~/components/plates/PlateCalculationExplanation.vue'
import type { Plate } from '~/types/lab'
import { countWellsWithoutLog10 } from '~/utils/plateDatasets'

/**
 * How a log10 measurement was calculated, e.g. for `label = 'Lum1_log10'`, `sourceLabel = 'Lum1'`.
 */
const props = defineProps<{
  plate: Plate
  label: string
  sourceLabel: string
  timestampIndex: number
}>()

const { t } = useI18n()

const texts = computed(() => {
  const skipped = countWellsWithoutLog10(props.plate, props.sourceLabel, props.label)
  return {
    title: t('plates.calculations.log10.title', { label: props.label }),
    formulas: [t('plates.calculations.log10.formula', { source: props.sourceLabel })],
    reading: t('plates.calculations.log10.reading'),
    warning: skipped > 0 ? t('plates.calculations.log10.skipped', { count: skipped }) : undefined,
  }
})
</script>

<template>
  <PlateCalculationExplanation
    :plate="props.plate"
    :label="props.label"
    :timestamp-index="props.timestampIndex"
    v-bind="texts"
  />
</template>
