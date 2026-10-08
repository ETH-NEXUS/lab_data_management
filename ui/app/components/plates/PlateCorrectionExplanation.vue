<script setup lang="ts">
import { computed } from 'vue'
import PlateDatasetSummary from '~/components/plates/PlateDatasetSummary.vue'
import type { Plate } from '~/types/lab'
import type { BackgroundCorrection } from '~/types/plateCalculations'
import { formatSummaryNumber, getWellTypeStatistic } from '~/utils/plateDatasets'

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

const methodName = computed(() => t(`plates.background_correction.methods.${props.correction.method}`))

// The value that was subtracted, e.g. the median of the R wells
const background = computed(() =>
  getWellTypeStatistic(
    props.plate,
    props.correction.source,
    props.correction.referenceType,
    props.correction.method,
    props.timestampIndex,
  ),
)
</script>

<template>
  <div class="mt-3">
    <h4 class="mb-1 font-medium text-slate-800">
      {{ t('plates.background_correction.heatmap_title', { label: props.correction.label }) }}
    </h4>
    <p class="text-sm text-slate-600">{{ t('plates.calculations.formula_caption') }}</p>
    <!-- Monospace with spaces, so "-" reads as a minus and not as a dash -->
    <code class="my-1 block w-fit rounded-md bg-slate-100 px-3 py-2 font-mono text-sm text-slate-800">
      {{ t('plates.background_correction.formula', { method: methodName, reference: props.correction.referenceType }) }}
    </code>
    <p class="text-xs text-slate-500">
      {{ t('plates.background_correction.formula_note', { reference: props.correction.referenceType }) }}
    </p>
    <p class="text-xs text-slate-600">
      {{
        t('plates.calculations.background_value', {
          method: methodName,
          reference: props.correction.referenceType,
          value: formatSummaryNumber(background),
        })
      }}
    </p>
    <PlateDatasetSummary :plate="props.plate" :label="props.correction.label" :timestamp-index="props.timestampIndex" />
  </div>
</template>
