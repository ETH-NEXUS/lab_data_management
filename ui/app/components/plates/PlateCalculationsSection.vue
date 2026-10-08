<script setup lang="ts">
import { computed, nextTick, ref, watch } from 'vue'
import ColorLegend from '~/components/plates/ColorLegend.vue'
import PlateActivityExplanation from '~/components/plates/PlateActivityExplanation.vue'
import PlateCalculationModal from '~/components/plates/PlateCalculationModal.vue'
import PlateCorrectionExplanation from '~/components/plates/PlateCorrectionExplanation.vue'
import PlateLog10Explanation from '~/components/plates/PlateLog10Explanation.vue'
import PlateTable from '~/components/plates/PlateTable.vue'
import { usePlateStore } from '~/stores/plates'
import { usePlateViewStore } from '~/stores/plateView'
import type { Plate, WellInfo } from '~/types/lab'
import type { PlateCalculationResult } from '~/types/plateCalculations'
import { getHeatmapRange } from '~/utils/heatmapScale'
import {
  findBackgroundCorrections,
  getActivityDataset,
  getCorrectionDataset,
  getLog10SourceLabel,
} from '~/utils/plateDatasets'

/**
 * The "Calculations" of the plate page: the button to start one, and below the
 * main heatmap how its measurement was calculated (log10, %Activity) or a
 * heatmap of its background corrections.
 */
const props = defineProps<{
  plate: Plate
}>()

const emit = defineEmits<{
  (e: 'well-selected', payload: WellInfo): void
}>()

const { t } = useI18n()
const plateStore = usePlateStore()
const plateViewStore = usePlateViewStore()

const isModalOpen = ref(false)

// The measurement of the main heatmap, if it is a log10, e.g. { label: 'Lum1_log10', source: 'Lum1' }
const log10Dataset = computed(() => {
  const label = plateViewStore.selectedMeasurement
  const source = getLog10SourceLabel(props.plate, label)
  return label && source ? { label, source } : null
})

// The measurement of the main heatmap, if it is a %Activity
const activityDataset = computed(() => getActivityDataset(props.plate, plateViewStore.selectedMeasurement))

// The measurement of the main heatmap, if it is a background correction itself
const correctionDataset = computed(() => getCorrectionDataset(props.plate, plateViewStore.selectedMeasurement))

// The background corrections of the measurement of the main heatmap
const corrections = computed(() => findBackgroundCorrections(props.plate, plateViewStore.selectedMeasurement))

const selectedLabel = ref<string | null>(null)
// A correction to show as soon as the reloaded plate has it, e.g. 'Lum1_bc_R_mean'
const labelToShow = ref<string | null>(null)

watch(
  corrections,
  () => {
    const labels = corrections.value.map((correction) => correction.label)
    if (labelToShow.value && labels.includes(labelToShow.value)) {
      selectedLabel.value = labelToShow.value
      labelToShow.value = null
    } else if (!selectedLabel.value || !labels.includes(selectedLabel.value)) {
      selectedLabel.value = labels[0] ?? null
    }
  },
  { immediate: true },
)

const selectedCorrection = computed(() => {
  return corrections.value.find((correction) => correction.label === selectedLabel.value) ?? null
})

// The corrected measurement has the same time points as the original one
const heatmapRange = computed(() =>
  getHeatmapRange(
    props.plate.details.overall_stats,
    selectedLabel.value,
    plateViewStore.selectedTimestampIdx,
    plateViewStore.heatmapScale,
  ),
)

/**
 * Shows a measurement in the main heatmap at the time point shown before:
 * DynamicPlate starts another measurement at its first time point, but a
 * calculated one has the same time points as the one it was calculated from.
 */
const showInMainHeatmap = async (label: string): Promise<void> => {
  const timestampIndex = plateViewStore.selectedTimestampIdx
  plateViewStore.selectedMeasurement = label
  await nextTick()
  plateViewStore.selectedTimestampIdx = timestampIndex
}

/**
 * Reloads only the plate (not the whole page), so the heatmap settings stay, and
 * shows the result: a log10 or %Activity in the main heatmap, a correction in
 * its own heatmap below (the main heatmap shows the measurement it is of).
 *
 * Accepted input example: `{ calculation: 'log10', label: 'Lum1', newLabel: 'Lum1_log10' }`
 */
const onCalculated = async (result: PlateCalculationResult): Promise<void> => {
  const isCorrection = result.calculation === 'background_correction'
  // Also when another correction of this measurement was shown before
  labelToShow.value = isCorrection ? result.newLabel : null
  try {
    const plate = await plateStore.fetchPlateByBarcode(props.plate.barcode)
    if (!plate) {
      labelToShow.value = null
      return
    }

    plateViewStore.measurementOptions = plate.details.measurement_labels ?? []
    // The well details on the right show the new measurement too
    const selectedPosition = plateViewStore.selectedWellInfo?.position
    if (selectedPosition !== undefined) {
      const selectedWell = plate.wells.find((well) => well.position === selectedPosition)
      plateViewStore.selectedWellInfo = { well: selectedWell, position: selectedPosition }
    }
    if (isCorrection) {
      await showInMainHeatmap(result.label)
    } else {
      plateViewStore.showHeatmap = true
      await showInMainHeatmap(result.newLabel)
    }
  } catch (err) {
    labelToShow.value = null
    // The plate store keeps the error and the page shows it
    console.error(err)
  }
}
</script>

<template>
  <section class="mt-8 rounded-xl border border-black/10 bg-white/50 p-4">
    <div class="flex flex-wrap items-center justify-between gap-3">
      <h3 class="text-lg font-medium text-slate-800">{{ t('plates.calculations.section_title') }}</h3>
      <UButton
        color="secondary"
        variant="outline"
        icon="i-heroicons-adjustments-horizontal"
        :label="t('plates.calculations.open_button')"
        @click="isModalOpen = true"
      />
    </div>

    <PlateLog10Explanation
      v-if="log10Dataset"
      :plate="props.plate"
      :label="log10Dataset.label"
      :source-label="log10Dataset.source"
      :timestamp-index="plateViewStore.selectedTimestampIdx"
    />

    <PlateActivityExplanation
      v-if="activityDataset"
      :plate="props.plate"
      :dataset="activityDataset"
      :timestamp-index="plateViewStore.selectedTimestampIdx"
    />

    <PlateCorrectionExplanation
      v-if="correctionDataset"
      :plate="props.plate"
      :correction="correctionDataset"
      :timestamp-index="plateViewStore.selectedTimestampIdx"
    />

    <p
      v-if="!log10Dataset && !activityDataset && !correctionDataset && !selectedCorrection"
      class="mt-2 text-sm text-slate-600"
    >
      {{ t('plates.calculations.none_yet', { label: plateViewStore.selectedMeasurement ?? '' }) }}
    </p>

    <template v-if="selectedCorrection">
      <PlateCorrectionExplanation
        class="mb-3"
        :plate="props.plate"
        :correction="selectedCorrection"
        :timestamp-index="plateViewStore.selectedTimestampIdx"
      />

      <div v-if="corrections.length > 1" class="mb-3 max-w-sm">
        <select
          v-model="selectedLabel"
          class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-lime-500"
        >
          <option v-for="correction in corrections" :key="`correction-${correction.label}`" :value="correction.label">
            {{ correction.label }}
          </option>
        </select>
      </div>

      <!-- Always a heatmap, also without "Show heatmap" for the main one -->
      <div class="flex flex-nowrap gap-4">
        <div class="min-w-0 overflow-auto">
          <PlateTable
            :plate="props.plate"
            :min="heatmapRange.min"
            :max="heatmapRange.max"
            :measurement-label="selectedCorrection.label"
            @well-selected="emit('well-selected', $event)"
          />
        </div>

        <ColorLegend :range="heatmapRange" always-shown />
      </div>
    </template>

    <PlateCalculationModal v-model:open="isModalOpen" :plate="props.plate" @calculated="onCalculated" />
  </section>
</template>
