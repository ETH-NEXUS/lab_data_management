<script setup lang="ts">
import { computed } from 'vue'
import HeatmapScaleNote from '~/components/plates/HeatmapScaleNote.vue'
import PlateStatsExplanation from '~/components/plates/PlateStatsExplanation.vue'
import { usePlateViewStore } from '~/stores/plateView'
import type { Plate } from '~/types/lab'
import type { HeatmapRange } from '~/utils/heatmapScale'
import { getCorrectionDataset, getLog10SourceLabel, getNormalizedDataset } from '~/utils/plateDatasets'

/**
 * One explanation above the heatmaps of all plates of an experiment (not one per
 * plate): which measurement and read they show, what its values are, whether
 * the colours can be compared between plates, and what z' and SSMD mean.
 *
 * Accepted props example:
 * - `{ plates: [...], timestamps: { Lum1: ['2026-09-30T15:48:08+02:00'] },
 *     range: { min: 33, max: 126340, lowerClipped: false, upperClipped: false } }`
 */
const props = defineProps<{
  plates: Plate[]
  timestamps: { [label: string]: string[] }
  // The colour scale of all plates together
  range: HeatmapRange
}>()

const { t } = useI18n()
const plateViewStore = usePlateViewStore()

const label = computed(() => plateViewStore.selectedMeasurement)

// The plates that have the measurement: one calculated on a plate page is only on that plate
const platesWithLabel = computed(() => {
  if (!label.value) return []
  const measurement = label.value
  return props.plates.filter((plate) => (plate.details.measurement_labels ?? []).includes(measurement))
})

// e.g. "Lum1, read 1 of 1 (2026-09-30T15:48:08+02:00)"
const title = computed(() => {
  if (!label.value) return ''
  const timestamps = props.timestamps[label.value] ?? []
  return t('experiments.results.heatmap.title', {
    label: label.value,
    number: plateViewStore.selectedTimestampIdx + 1,
    total: timestamps.length,
    time: timestamps[plateViewStore.selectedTimestampIdx] ?? '–',
  })
})

/**
 * What the values of the measurement are, in one sentence. A calculation of a
 * plate page is recognized by its name on a plate that has it, as on the plate page.
 */
const meaning = computed(() => {
  const plate = platesWithLabel.value[0]
  if (!label.value || !plate) return ''

  const log10Source = getLog10SourceLabel(plate, label.value)
  if (log10Source) {
    return t('experiments.results.heatmap.meaning.log10', { source: log10Source })
  }
  const correction = getCorrectionDataset(plate, label.value)
  if (correction) {
    return t('experiments.results.heatmap.meaning.correction', {
      source: correction.source,
      reference: correction.referenceType,
      method: t(`plates.background_correction.methods.${correction.method}`),
    })
  }
  const normalized = getNormalizedDataset(plate, label.value)
  if (normalized) {
    return t(`experiments.results.heatmap.meaning.${normalized.kind}`, {
      source: normalized.source,
      negative: normalized.negativeType,
      positive: normalized.positiveType,
    })
  }
  return t('experiments.results.heatmap.meaning.raw')
})

const missingPlates = computed(() => props.plates.length - platesWithLabel.value.length)
</script>

<template>
  <div v-if="label" class="mb-6 max-w-4xl space-y-1 rounded-xl border border-black/10 bg-white/50 p-4">
    <h4 class="font-medium text-slate-800">{{ title }}</h4>
    <p class="text-sm text-slate-700">{{ meaning }}</p>
    <p class="text-sm text-slate-700">{{ t('experiments.results.heatmap.open_plate') }}</p>
    <p v-if="missingPlates > 0" class="text-xs text-amber-700">
      {{ t('experiments.results.heatmap.missing_plates', { count: missingPlates, total: props.plates.length }) }}
    </p>

    <template v-if="plateViewStore.perPlateView">
      <p class="text-sm text-slate-700">{{ t('experiments.results.heatmap.scale_per_plate') }}</p>
    </template>
    <template v-else>
      <p class="text-sm text-slate-700">{{ t('experiments.results.heatmap.scale_shared') }}</p>
      <HeatmapScaleNote :range="props.range" />
    </template>

    <PlateStatsExplanation
      v-if="plateViewStore.selectedPosControl && plateViewStore.selectedNegControl"
      :positive="plateViewStore.selectedPosControl"
      :negative="plateViewStore.selectedNegControl"
    />
  </div>
</template>
