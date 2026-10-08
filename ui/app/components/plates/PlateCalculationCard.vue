<script setup lang="ts">
import { computed } from 'vue'
import ColorLegend from '~/components/plates/ColorLegend.vue'
import HeatmapScaleNote from '~/components/plates/HeatmapScaleNote.vue'
import PlateCalculationDeleteButton from '~/components/plates/PlateCalculationDeleteButton.vue'
import PlateCorrectionExplanation from '~/components/plates/PlateCorrectionExplanation.vue'
import PlateLog10Explanation from '~/components/plates/PlateLog10Explanation.vue'
import PlateNormalizationExplanation from '~/components/plates/PlateNormalizationExplanation.vue'
import PlateTable from '~/components/plates/PlateTable.vue'
import { usePlateViewStore } from '~/stores/plateView'
import type { Plate, WellInfo } from '~/types/lab'
import { getHeatmapRange } from '~/utils/heatmapScale'
import { getCorrectionDataset, getLog10SourceLabel, getNormalizedDataset } from '~/utils/plateDatasets'

/**
 * One measurement calculated on the plate page, e.g. `Lum1_bc_R_median`: how it
 * was calculated, its own heatmap (always shown, also without "Show heatmap"),
 * and a button to delete it.
 */
const props = defineProps<{
  plate: Plate
  label: string
}>()

const emit = defineEmits<{
  (e: 'well-selected', payload: WellInfo): void
  (e: 'deleted', label: string): void
}>()

const plateViewStore = usePlateViewStore()

// The read chosen above for the main heatmap; the explanation names it
const timestampIndex = computed(() => plateViewStore.selectedTimestampIdx)

// How the measurement was calculated: one of these is set
const log10Source = computed(() => getLog10SourceLabel(props.plate, props.label))
const normalizedDataset = computed(() => getNormalizedDataset(props.plate, props.label))
const correctionDataset = computed(() => getCorrectionDataset(props.plate, props.label))

const heatmapRange = computed(() =>
  getHeatmapRange(props.plate.details.overall_stats, props.label, timestampIndex.value, plateViewStore.heatmapScale),
)
</script>

<template>
  <article class="mt-6 rounded-lg border border-black/10 bg-white/60 p-4">
    <div class="flex justify-end">
      <PlateCalculationDeleteButton :plate="props.plate" :label="props.label" @deleted="emit('deleted', $event)" />
    </div>

    <PlateLog10Explanation
      v-if="log10Source"
      :plate="props.plate"
      :label="props.label"
      :source-label="log10Source"
      :timestamp-index="timestampIndex"
    />
    <PlateNormalizationExplanation
      v-else-if="normalizedDataset"
      :plate="props.plate"
      :dataset="normalizedDataset"
      :timestamp-index="timestampIndex"
    />
    <PlateCorrectionExplanation
      v-else-if="correctionDataset"
      :plate="props.plate"
      :correction="correctionDataset"
      :timestamp-index="timestampIndex"
    />

    <div class="mt-3 flex flex-nowrap gap-4">
      <div class="min-w-0 overflow-auto">
        <PlateTable
          :plate="props.plate"
          :min="heatmapRange.min"
          :max="heatmapRange.max"
          :measurement-label="props.label"
          @well-selected="emit('well-selected', $event)"
        />
      </div>
      <ColorLegend :range="heatmapRange" always-shown />
    </div>
    <HeatmapScaleNote :range="heatmapRange" />
  </article>
</template>
