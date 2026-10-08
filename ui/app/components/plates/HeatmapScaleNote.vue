<script setup lang="ts">
import { computed } from 'vue'
import { usePlateViewStore } from '~/stores/plateView'
import type { HeatmapRange } from '~/utils/heatmapScale'
import { formatSummaryNumber } from '~/utils/plateDatasets'

/**
 * One line below a heatmap: how the values become colours, with the values at
 * the two ends of the scale, e.g. "Colours: from the lowest value (33) ...".
 */
const props = defineProps<{
  range: HeatmapRange
}>()

const { t } = useI18n()
const plateViewStore = usePlateViewStore()

const text = computed(() =>
  t(`plates.calculations.scale_note.${plateViewStore.heatmapScale}`, {
    min: formatSummaryNumber(props.range.min),
    max: formatSummaryNumber(props.range.max),
  }),
)
</script>

<template>
  <p class="mt-1 max-w-3xl text-xs text-slate-500">{{ text }}</p>
</template>
