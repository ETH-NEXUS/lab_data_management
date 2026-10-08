<script setup lang="ts">
import { computed } from 'vue'
import { usePlateViewStore } from '~/stores/plateView'
import type { HeatmapRange } from '~/utils/heatmapScale'
import { buildPlateLegend } from '~/utils/plateHeatmap'

type Props = {
  // e.g. { min: 33, max: 1934, lowerClipped: false, upperClipped: true }
  range: HeatmapRange
  // Shown also without "Show heatmap", e.g. next to the background corrected heatmap
  alwaysShown?: boolean
}

const props = defineProps<Props>()

// The legend has 21 colours (0..20), every fifth of them with its value
const LEGEND_STEPS = 20
const LABEL_EVERY = 5
const platePage = usePlateViewStore()

const legendColors = computed(() => {
  if (!platePage.selectedMeasurement) return undefined
  return buildPlateLegend(props.range.min, props.range.max, platePage.heatmapPalette, LEGEND_STEPS)
})

/**
 * The text next to a step of the legend; only every fifth step has one. The
 * legend runs from the top (max) to the bottom (min), and a clipped end also
 * stands for the wells beyond it: "≥ 1934.0", "≤ 33.0".
 */
const stepLabel = (value: number, index: number): string => {
  if (index % LABEL_EVERY !== 0) return ' '
  const text = value.toFixed(1)
  if (index === 0 && props.range.upperClipped) return `≥ ${text}`
  if (index === LEGEND_STEPS && props.range.lowerClipped) return `≤ ${text}`
  return text
}
</script>

<template>
  <div
    v-if="(props.alwaysShown || platePage.showHeatmap) && platePage.selectedMeasurement && legendColors"
    class="legendWrap"
  >
    <div
      v-for="(color, idx) in legendColors"
      :key="color.value + idx"
      class="legendItem"
      :style="{ backgroundColor: color.color }"
    >
      <span class="legendLabel">{{ stepLabel(color.value, idx) }}</span>
    </div>
  </div>
</template>

<style scoped>
.legendWrap {
  margin-top: 1rem;
  margin-bottom: 1rem;
  margin-left: 1rem;
}

.legendItem {
  position: relative;
  width: 30px;
  height: 10px;
}

.legendLabel {
  position: absolute;
  left: 33px;
  font-size: 9px;
  white-space: nowrap;
}
</style>
