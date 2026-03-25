<script setup lang="ts">
import { computed } from 'vue'
import type { Well, WellDetails } from '~/types/lab'
import { usePlateViewStore } from '~/stores/plateView'
const { t } = useI18n()

const props = defineProps<{
  well: Well | WellDetails | undefined
  row: number
  col: number
}>()

const plateViewStore = usePlateViewStore()

const measurement = (well: WellDetails) => {
  const selectedMeasurement = plateViewStore.selectedMeasurement
  const selectedTimestampIdx = plateViewStore.selectedTimestampIdx

  if (!selectedMeasurement || !well.measurements) {
    return null
  }

  if (!(selectedMeasurement in well.measurements)) {
    return null
  }

  const series = well.measurements[selectedMeasurement]
  if (!series || selectedTimestampIdx < 0 || series.length <= selectedTimestampIdx) {
    return null
  }

  return series[selectedTimestampIdx] ?? null
}

const compoundLabels = computed(() => {
  if (!props.well || !('compounds' in props.well) || !Array.isArray(props.well.compounds)) {
    return [] as string[]
  }

  return props.well.compounds
    .map((compound) => {
      if (typeof compound === 'string') return compound
      if (compound && typeof compound === 'object' && 'name' in compound && typeof compound.name === 'string') {
        return compound.name
      }
      return null
    })
    .filter((label): label is string => Boolean(label))
})

const measurementValue = computed(() => {
  if (!props.well || !('measurements' in props.well)) return null
  return measurement(props.well as WellDetails)
})

const tooltipText = computed(() => {
  if (!props.well) {
    return ''
  }

  const lines: string[] = []
  lines.push(`${props.well.hr_position} (${props.well.type}) (${props.well.position})`)

  if (compoundLabels.value.length > 0) {
    lines.push(t('label.compounds'))
    lines.push(...compoundLabels.value)
  }

  if (plateViewStore.selectedMeasurement) {
    lines.push(t('label.measurements'))
    lines.push(`${plateViewStore.selectedMeasurement}: ${measurementValue.value ?? 'N/A'}`)
  }

  return lines.join('\n')
})
</script>

<template>
  <span v-if="props.well" class="inline-flex h-full w-full items-center justify-center" :title="tooltipText">
    <slot />
  </span>
  <slot v-else />
</template>
