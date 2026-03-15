<script setup lang="ts">
import { computed } from 'vue'

const props = defineProps<{
  name?: string
  series?: Array<number | string | Date>
  categories?: Array<string | Date | number>
  color?: string
}>()

const chartWidth = 640
const chartHeight = 250
const padding = { top: 20, right: 20, bottom: 30, left: 40 }

const numericSeries = computed<number[]>(() =>
  (props.series ?? [])
    .map((value) => (typeof value === 'string' ? Number.parseFloat(value) : value))
    .filter((value): value is number => Number.isFinite(value)),
)

const hasChart = computed(() => {
  const categories = props.categories ?? []
  return categories.length > 1 && numericSeries.value.length > 1
})

const minValue = computed(() => Math.min(...numericSeries.value))
const maxValue = computed(() => Math.max(...numericSeries.value))

const points = computed(() => {
  if (!hasChart.value) return ''

  const series = numericSeries.value
  const xStep = (chartWidth - padding.left - padding.right) / Math.max(series.length - 1, 1)
  const valueRange = maxValue.value - minValue.value || 1

  return series
    .map((value, index) => {
      const x = padding.left + index * xStep
      const normalized = (value - minValue.value) / valueRange
      const y = chartHeight - padding.bottom - normalized * (chartHeight - padding.top - padding.bottom)
      return `${x},${y}`
    })
    .join(' ')
})

const strokeColor = computed(() => props.color || '#e08214')
</script>

<template>
  <div :id="`chart-${props.name}`" class="w-full">
    <p class="mb-2 text-sm font-semibold text-slate-700">{{ props.name }}</p>
    <svg v-if="hasChart" :viewBox="`0 0 ${chartWidth} ${chartHeight}`" class="h-[250px] w-full">
      <line
        :x1="padding.left"
        :y1="chartHeight - padding.bottom"
        :x2="chartWidth - padding.right"
        :y2="chartHeight - padding.bottom"
        stroke="#e7e7e7"
        stroke-width="1"
      />
      <line
        :x1="padding.left"
        :y1="padding.top"
        :x2="padding.left"
        :y2="chartHeight - padding.bottom"
        stroke="#e7e7e7"
        stroke-width="1"
      />
      <polyline fill="none" :stroke="strokeColor" stroke-width="2" :points="points" />
    </svg>
  </div>
</template>
