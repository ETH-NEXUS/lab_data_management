<script setup lang="ts">
import { computed } from 'vue'
import { usePlateViewStore } from '~/stores/plateView'
import { buildPlateLegend } from '~/utils/plateHeatmap'

type Props = {
  min: number
  max: number
}

const props = defineProps<Props>()
const platePage = usePlateViewStore()

const legendColors = computed(() => {
  if (!platePage.selectedMeasurement) return undefined
  return buildPlateLegend(props.min, props.max, platePage.heatmapPalette, 20)
})
</script>

<template>
  <div v-if="platePage.showHeatmap && platePage.selectedMeasurement && legendColors" class="legendWrap">
    <div
      v-for="(color, idx) in legendColors"
      :key="color.value + idx"
      class="legendItem"
      :style="{ backgroundColor: color.color }"
    >
      <span class="legendLabel">{{ [0, 5, 10, 15, 20].includes(idx) ? color.value.toFixed(1) : ' ' }}</span>
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
