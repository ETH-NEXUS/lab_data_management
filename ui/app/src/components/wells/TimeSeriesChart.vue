<script setup lang="ts">
import {defineProps, onMounted, reactive} from 'vue'
import ApexCharts from 'apexcharts'
import {TimeSeriesChart} from 'components/models'

const props = defineProps({
  name: {
    type: String,
  },
  series: {
    type: Array<number>,
  },
  categories: {
    type: Array<string | Date>,
  },
  color: {
    type: String,
  },
})

onMounted(async () => {
  createOptions()
  if (options.xaxis.categories && options.xaxis.categories.length > 1) {
    const chart = new ApexCharts(document.querySelector(`#chart-${props.name}`), options)
    await chart.render()
  }
})

const options = reactive<TimeSeriesChart>({
  series: [],
  colors: ['#e08214'],
  chart: {
    id: 'time-series',
    type: 'line',
    height: 250,
    zoom: {
      enabled: false,
    },
    dropShadow: {
      enabled: true,
      color: '#000',
      top: 18,
      left: 7,
      blur: 10,
      opacity: 0.2,
    },
  },
  dataLabels: {
    enabled: false,
  },
  stroke: {
    curve: 'straight',
  },
  title: {
    text: 'Time Series',
    align: 'left',
  },
  grid: {
    borderColor: '#e7e7e7',
    row: {
      colors: ['#f3f3f3', 'transparent'], // takes an array which will be repeated on columns
      opacity: 0.5,
    },
  },
  xaxis: {
    categories: [],
    type: 'datetime',
  },
  yaxis: {
    title: {
      text: 'Value',
    },
  },
  markers: {
    size: 1,
  },

  legend: {
    position: 'top',
    horizontalAlign: 'left',
    floating: true,
    offsetY: -25,
    offsetX: -5,
  },
})

const createOptions = () => {
  options.series.push({
    name: props.name,
    data: props.series,
  })
  options.xaxis.categories = props.categories
  options.title.text = props.name

  if (props.color) {
    options.colors = [props.color]
  }
}
</script>

<template>
  <div :id="`chart-${props.name}`"></div>
</template>

<style lang="sass"></style>
