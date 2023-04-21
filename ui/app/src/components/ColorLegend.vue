<script setup lang="ts">
import {storeToRefs} from 'pinia'
import {useSettingsStore} from 'stores/settings'
import {computed, defineProps} from 'vue'
import {LegendColor} from 'components/models'
import {percentageToHsl} from 'components/helpers'

const props = defineProps({
  min: {
    type: Number,
    required: true,
  },
  max: {
    type: Number,
    required: true,
  },
  selectedMeasurement: {
    default: null,
  },
})

const legendColors = computed(() => {
  return createColorLegend()
})

const {platePage} = storeToRefs(useSettingsStore())

const createColorLegend = () => {
  const legend: LegendColor[] = []
  if (props.selectedMeasurement) {
    const numberOfSteps = 10

    const step: number = (props.max - props.min) / numberOfSteps
    for (let i = 0; i <= numberOfSteps; i++) {
      const value: number = props.min + i * step
      const color = percentageToHsl(
        (value - props.min) / (props.max - props.min),
        platePage.value.heatmapPalette.value.from,
        platePage.value.heatmapPalette.value.to
      )
      legend.push({value, color})
    }

    return legend.reverse()
  }
}
</script>

<template>
  <div v-if="platePage.showHeatmap && selectedMeasurement && legendColors" class="q-my-md q-ml-md">
    <div
      class="legendItem"
      v-for="(color, idx) in legendColors"
      :key="color.value + idx"
      :style="{backgroundColor: color.color}">
      <span class="legendLabel">{{ [0, 5, 10].includes(idx) ? color.value.toFixed(3) : ' ' }}</span>
    </div>
  </div>
</template>

<style scoped lang="sass">
.legendItem
  position: relative
  width: 30px
  height: 15px

.legendLabel
  position: absolute
  left: 33px
  font-size: 10px
</style>
