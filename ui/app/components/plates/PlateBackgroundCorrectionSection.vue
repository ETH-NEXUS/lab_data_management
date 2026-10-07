<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import ColorLegend from '~/components/plates/ColorLegend.vue'
import PlateBackgroundCorrectionModal from '~/components/plates/PlateBackgroundCorrectionModal.vue'
import PlateActivityExplanation from '~/components/plates/PlateActivityExplanation.vue'
import PlateDatasetSummary from '~/components/plates/PlateDatasetSummary.vue'
import PlateTable from '~/components/plates/PlateTable.vue'
import { usePlateStore } from '~/stores/plates'
import { usePlateViewStore } from '~/stores/plateView'
import type { PlateCalculationResult } from '~/types/backgroundCorrection'
import type { Plate, WellInfo } from '~/types/lab'
import { findBackgroundCorrections } from '~/utils/backgroundCorrection'
import {
  countWellsWithoutLog10,
  formatSummaryNumber,
  getActivityDataset,
  getLog10SourceLabel,
} from '~/utils/plateDatasets'
import { getHeatmapRange } from '~/utils/heatmapScale'

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

// The corrected versions of the measurement chosen for the main heatmap
const corrections = computed(() => findBackgroundCorrections(props.plate, plateViewStore.selectedMeasurement))

const selectedLabel = ref<string | null>(null)
// A correction to show as soon as the reloaded plate has it, e.g. 'Lum_CTG_bc_N1_mean'
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

// The measurement chosen for the main heatmap, if it is the log10 of another one, e.g. 'Lum1'
const log10Source = computed(() => getLog10SourceLabel(props.plate, plateViewStore.selectedMeasurement))

// The measurement chosen for the main heatmap, if it is a %Activity
const activityDataset = computed(() => getActivityDataset(props.plate, plateViewStore.selectedMeasurement))

const wellsWithoutLog10 = computed(() => {
  if (!log10Source.value || !plateViewStore.selectedMeasurement) return 0
  return countWellsWithoutLog10(props.plate, log10Source.value, plateViewStore.selectedMeasurement)
})

// The value that was subtracted, e.g. the median of the R wells at the shown time point
const background = computed(() => {
  const correction = selectedCorrection.value
  const label = plateViewStore.selectedMeasurement
  if (!correction || !label) return null
  const statsOfReference = props.plate.details.stats[label]?.[correction.referenceType]
  return statsOfReference?.[correction.method]?.[plateViewStore.selectedTimestampIdx] ?? null
})

// The corrected measurement has the same time points as the original one
const heatmapRange = computed(() => {
  const stats = selectedLabel.value ? props.plate.details.overall_stats[selectedLabel.value] : undefined
  return getHeatmapRange(stats, plateViewStore.selectedTimestampIdx, plateViewStore.heatmapScale)
})

/**
 * Reloads only the plate (not the whole page), so the heatmap settings stay, and
 * shows the result: a log10 in the main heatmap, a correction below it (the main
 * heatmap then shows the measurement it was calculated from).
 *
 * Accepted input example: `{ calculation: 'log10', label: 'Lum1', newLabel: 'Lum1_log10' }`
 */
const onCalculated = async (result: PlateCalculationResult): Promise<void> => {
  if (result.calculation === 'background_correction') {
    // Also when another correction of this measurement was shown before
    labelToShow.value = result.newLabel
  }
  try {
    const plate = await plateStore.fetchPlateByBarcode(props.plate.barcode)
    if (!plate) return

    plateViewStore.measurementOptions = plate.details.measurement_labels ?? []
    // A log10 or a %Activity replaces the measurement of the main heatmap
    if (result.calculation === 'log10' || result.calculation === 'percent_activity') {
      plateViewStore.selectedMeasurement = result.newLabel
      plateViewStore.showHeatmap = true
    } else {
      plateViewStore.selectedMeasurement = result.label
    }
  } catch (err) {
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

    <!-- The main heatmap shows a log10: how it was calculated and what came out -->
    <div v-if="log10Source && plateViewStore.selectedMeasurement" class="mt-3">
      <h4 class="mb-1 font-medium text-slate-800">
        {{ t('plates.calculations.log10_title', { label: plateViewStore.selectedMeasurement }) }}
      </h4>
      <p class="text-sm text-slate-600">{{ t('plates.background_correction.formula_caption') }}</p>
      <code class="my-1 block w-fit rounded-md bg-slate-100 px-3 py-2 font-mono text-sm text-slate-800">
        {{ t('plates.calculations.log10_formula', { label: log10Source }) }}
      </code>
      <p v-if="wellsWithoutLog10 > 0" class="text-xs text-amber-700">
        {{ t('plates.calculations.log10_skipped', { count: wellsWithoutLog10 }) }}
      </p>
      <PlateDatasetSummary
        :plate="props.plate"
        :label="plateViewStore.selectedMeasurement"
        :timestamp-index="plateViewStore.selectedTimestampIdx"
      />
    </div>

    <!-- The main heatmap shows a %Activity: its formula with the numbers put in -->
    <div v-if="activityDataset" class="mt-3">
      <h4 class="mb-1 font-medium text-slate-800">
        {{ t('plates.calculations.activity_title', { label: activityDataset.label }) }}
      </h4>
      <p class="text-sm text-slate-600">{{ t('plates.background_correction.formula_caption') }}</p>
      <PlateActivityExplanation
        :plate="props.plate"
        :dataset="activityDataset"
        :timestamp-index="plateViewStore.selectedTimestampIdx"
      />
      <PlateDatasetSummary
        :plate="props.plate"
        :label="activityDataset.label"
        :timestamp-index="plateViewStore.selectedTimestampIdx"
      />
    </div>

    <p v-if="!selectedCorrection && !log10Source && !activityDataset" class="mt-2 text-sm text-slate-600">
      {{ t('plates.calculations.none_yet', { label: plateViewStore.selectedMeasurement ?? '' }) }}
    </p>

    <template v-if="selectedCorrection">
      <h4 class="mt-3 mb-1 font-medium text-slate-800">
        {{ t('plates.background_correction.heatmap_title', { label: selectedCorrection.label }) }}
      </h4>
      <p class="text-sm text-slate-600">{{ t('plates.background_correction.formula_caption') }}</p>
      <!-- Monospace with spaces, so "-" reads as a minus and not as a dash -->
      <code class="my-1 block w-fit rounded-md bg-slate-100 px-3 py-2 font-mono text-sm text-slate-800">
        {{
          t('plates.background_correction.formula', {
            method: t(`plates.background_correction.methods.${selectedCorrection.method}`),
            reference: selectedCorrection.referenceType,
          })
        }}
      </code>
      <p class="text-xs text-slate-500">
        {{ t('plates.background_correction.formula_note', { reference: selectedCorrection.referenceType }) }}
      </p>
      <p class="text-xs text-slate-600">
        {{
          t('plates.calculations.background_value', {
            method: t(`plates.background_correction.methods.${selectedCorrection.method}`),
            reference: selectedCorrection.referenceType,
            value: formatSummaryNumber(background),
          })
        }}
      </p>
      <PlateDatasetSummary
        class="mb-3"
        :plate="props.plate"
        :label="selectedCorrection.label"
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

    <PlateBackgroundCorrectionModal v-model:open="isModalOpen" :plate="props.plate" @calculated="onCalculated" />
  </section>
</template>
