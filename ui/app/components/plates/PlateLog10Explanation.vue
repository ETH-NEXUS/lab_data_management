<script setup lang="ts">
import { computed } from 'vue'
import PlateDatasetSummary from '~/components/plates/PlateDatasetSummary.vue'
import type { Plate } from '~/types/lab'
import { countWellsWithoutLog10 } from '~/utils/plateDatasets'

/**
 * How a log10 measurement was calculated and what came out, e.g. for
 * `label = 'Lum1_log10'`, `sourceLabel = 'Lum1'`.
 */
const props = defineProps<{
  plate: Plate
  label: string
  sourceLabel: string
  timestampIndex: number
}>()

const { t } = useI18n()

const wellsWithoutLog10 = computed(() => countWellsWithoutLog10(props.plate, props.sourceLabel, props.label))
</script>

<template>
  <div class="mt-3">
    <h4 class="mb-1 font-medium text-slate-800">{{ t('plates.calculations.log10_title', { label: props.label }) }}</h4>
    <p class="text-sm text-slate-600">{{ t('plates.calculations.formula_caption') }}</p>
    <code class="my-1 block w-fit rounded-md bg-slate-100 px-3 py-2 font-mono text-sm text-slate-800">
      {{ t('plates.calculations.log10_formula', { label: props.sourceLabel }) }}
    </code>
    <p v-if="wellsWithoutLog10 > 0" class="text-xs text-amber-700">
      {{ t('plates.calculations.log10_skipped', { count: wellsWithoutLog10 }) }}
    </p>
    <PlateDatasetSummary :plate="props.plate" :label="props.label" :timestamp-index="props.timestampIndex" />
  </div>
</template>
