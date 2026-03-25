<script setup lang="ts">
import { onMounted, ref } from 'vue'
import TimeSeriesChart from '~/components/wells/TimeSeriesChart.vue'
import type { Measurement, PlotData } from '~/types/lab'

const props = defineProps<{
  measurements?: Measurement[]
}>()

const predefinedColors = [
  '#67001f',
  '#4393c3',
  '#b2182b',
  '#f4a582',
  '#8c510a',
  '#7fbc41',
  '#de77ae',
  '#e08214',
  '#542788',
  '#878787',
  '#e6f598',
]

const plotData = ref<PlotData[]>([])

onMounted(async () => {
  createPlotData()
})

const createPlotData = () => {
  const dataByLabel = new Map<string, { x: string | Date; y: number }[]>()
  if (!props.measurements) {
    return
  }
  for (const measurement of props.measurements) {
    const label = measurement.label
    if (!dataByLabel.has(label)) {
      dataByLabel.set(label, [])
    }
    const data = dataByLabel.get(label)
    if (data) {
      data.push({ x: measurement.measured_at, y: measurement.value })
    }
  }
  const generatedPlotData: PlotData[] = []

  for (const [label, data] of dataByLabel) {
    generatedPlotData.push({
      name: label,
      data: data.map((entry) => entry.y.toFixed(2)),
      color: predefinedColors[generatedPlotData.length] ?? '#e08214',
      categories: data.map((entry) => entry.x),
    })
  }

  plotData.value = generatedPlotData
}
</script>

<template>
  <div v-for="item in plotData" :key="item.name" class="mt-4">
    <TimeSeriesChart :name="item.name" :series="item.data" :categories="item.categories" :color="item.color" />
  </div>
</template>
