<script setup lang="ts">
import {defineProps, onMounted, ref} from 'vue'
import {Measurement} from 'components/models'
import TimeSeriesChart from 'components/wells/TimeSeriesChart.vue'
import {PlotData} from 'components/models'

const props = defineProps({
  measurements: {
    type: Array<Measurement>,
  },
})
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

onMounted(async () => {
  createPlotData()
})
const plotData = ref<Array<PlotData>>([])

const createPlotData = () => {
  const dataByLabel = new Map<string, {x: string | Date; y: number}[]>()
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
      data.push({x: measurement.measured_at, y: measurement.value})
    }
  }
  const _plotData = []

  for (const [label, data] of dataByLabel) {
    _plotData.push({
      name: label,
      data: data.map(d => d.y.toFixed(2)),
      color: predefinedColors[_plotData.length],
      categories: data.map(d => d.x),
    })
  }

  plotData.value = _plotData
  console.log(JSON.stringify(plotData.value))
}
</script>

<template>
  <div v-for="item in plotData" :key="item.name" class="q-mt-lg">
    <TimeSeriesChart
      :name="item.name"
      :series="item.data"
      :categories="item.categories"
      :color="item.color" />
  </div>
</template>

<style lang="sass"></style>
