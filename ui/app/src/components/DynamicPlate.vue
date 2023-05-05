<script setup lang="ts">
import {computed, defineProps, defineEmits, PropType, ref, onMounted} from 'vue'
import {Plate, WellInfo} from './models'

import {palettes} from 'components/helpers'
import {storeToRefs} from 'pinia'
import {useSettingsStore} from 'stores/settings'
import {useI18n} from 'vue-i18n'
import MeasurementCalculator from 'components/MeasurementCalculator.vue'
import {useProjectStore} from 'stores/project'
import bus from 'src/eventBus'
import {useQuasar} from 'quasar'

import ColorLegend from 'components/ColorLegend.vue'
import HeatMapSettings from 'components/HeatMapSettings.vue'
import PlateStats from 'components/PlateStats.vue'
import PlateTable from 'components/PlateTable.vue'

const {t} = useI18n()

const {platePage} = storeToRefs(useSettingsStore())
const projectStore = useProjectStore()

const props = defineProps({
  plate: {
    type: Object as PropType<Plate>,
    required: true,
  },
})

const combinedLabels = computed<string[]>(() => {
  return props.plate.details.measurement_labels
})

const wellContentOptions = [
  {
    label: t('label.hr_position'),
    value: 'hr_position',
  },
  {
    label: t('label.well_type'),
    value: 'type',
  },
  {
    label: t('label.index_position'),
    value: 'position',
  },
]

const $q = useQuasar()

const selectedMeasurement = ref<string | null>(null)
const measurementOptions = ref<string[] | null>(null)

const selectedTimestampIdx = ref<number>(0)
const timestampOptions = computed(() => {
  if (selectedMeasurement.value) {
    const timestamps = props.plate.details.measurement_timestamps[selectedMeasurement.value]
    if (timestamps && timestamps.length > 0) {
      return timestamps.map((ts, idx) => {
        return {label: ts, value: idx}
      })
    }
  }
  return null
})

const openCalculator = ref<boolean>(false)

const emit = defineEmits(['well-selected', 'refresh'])

const z_prime = computed(() => {
  if (selectedMeasurement.value) {
    if (
      'P' in props.plate.details.stats[selectedMeasurement.value] &&
      'N' in props.plate.details.stats[selectedMeasurement.value]
    ) {
      const mad_pos =
        props.plate.details.stats[selectedMeasurement.value]['P'].mad[selectedTimestampIdx.value]
      const mad_neg =
        props.plate.details.stats[selectedMeasurement.value]['N'].mad[selectedTimestampIdx.value]
      const median_pos =
        props.plate.details.stats[selectedMeasurement.value]['P'].median[selectedTimestampIdx.value]
      const median_neg =
        props.plate.details.stats[selectedMeasurement.value]['N'].median[selectedTimestampIdx.value]
      return 1 - (3 * (mad_pos + mad_neg)) / Math.abs(median_pos - median_neg)
    }
  }
  return null
})

const ssmd = computed(() => {
  if (selectedMeasurement.value) {
    if (
      'P' in props.plate.details.stats[selectedMeasurement.value] &&
      'N' in props.plate.details.stats[selectedMeasurement.value]
    ) {
      const mad_pos =
        props.plate.details.stats[selectedMeasurement.value]['P'].mad[selectedTimestampIdx.value]
      const mad_neg =
        props.plate.details.stats[selectedMeasurement.value]['N'].mad[selectedTimestampIdx.value]
      const median_pos =
        props.plate.details.stats[selectedMeasurement.value]['P'].median[selectedTimestampIdx.value]
      const median_neg =
        props.plate.details.stats[selectedMeasurement.value]['N'].median[selectedTimestampIdx.value]
      return Math.abs(median_pos - median_neg) / (0.5 * (mad_pos + mad_neg))
    }
  }
  return null
})

const min = computed<number>(() => {
  if (!selectedMeasurement.value) return 0
  return props.plate.details.overall_stats[selectedMeasurement.value].min[selectedTimestampIdx.value]
})

const max = computed<number>(() => {
  if (!selectedMeasurement.value) return 0
  return props.plate.details.overall_stats[selectedMeasurement.value].max[selectedTimestampIdx.value]
})

onMounted(() => {
  if (props.plate.details.measurement_labels.length > 0) {
    measurementOptions.value = props.plate.details.measurement_labels
    selectedMeasurement.value = measurementOptions.value[0]
  }
  bus.on('openCalculator', () => {
    openCalculator.value = true
  })
})

const calculateNewMeasurement = async (expression: string, newLabel: string, usedLabels: string[]) => {
  openCalculator.value = false
  $q.loading.show({
    message: t('info.calculation_in_progress'),
  })
  await projectStore.addNewMeasurement(props.plate.id, expression, newLabel, usedLabels)
  $q.loading.hide()

  emit('refresh')
}
</script>

<template>
  <HeatMapSettings />
  <div class="row no-wrap">
    <PlateTable
      @well-selected="(well_info: WellInfo) => (platePage.selectedWellInfo = well_info)"
      :max="max"
      :min="min"
      :plate="plate"
      :selected-measurement="selectedMeasurement"
      :selected-timestamp-idx="selectedTimestampIdx" />

    <ColorLegend :max="max" :min="min" :selectedMeasurement="selectedMeasurement" />
  </div>
  <div class="row">
    <PlateStats :ssmd="ssmd" :z_prime="z_prime" v-if="ssmd && z_prime" />
  </div>

  <div v-if="!plate.template" class="row">
    <div class="col-4 q-mr-md">
      <q-select
        v-model="platePage.wellContent"
        :options="wellContentOptions"
        :label="t('label.well_content')"
        emit-value
        map-options />
    </div>

    <div v-if="measurementOptions" class="col-4">
      <q-checkbox
        v-model="platePage.showHeatmap"
        :label="t('label.show_heatmap')"
        class="q-mr-md"></q-checkbox>

      <div class="col-12" v-if="platePage.showHeatmap && measurementOptions">
        <q-select
          v-if="measurementOptions.length > 0"
          v-model="selectedMeasurement"
          :disable="measurementOptions.length === 1"
          :options="measurementOptions"
          :label="t('message.select_label')"></q-select>
        <q-select
          v-if="timestampOptions"
          v-model="selectedTimestampIdx"
          :disable="timestampOptions.length === 1"
          :options="timestampOptions"
          :label="t('message.select_timestamp')"
          emit-value
          map-options></q-select>
      </div>

      <div class="col-12" v-if="platePage.showHeatmap && measurementOptions.length > 0">
        <q-select
          v-model="platePage.heatmapPalette"
          :options="palettes"
          :label="t('label.select_color_palette')"></q-select>
      </div>
    </div>
  </div>
  <q-dialog v-model="openCalculator">
    <q-card class="calculator_dialog">
      <q-card-section>
        <MeasurementCalculator
          :labels="combinedLabels"
          @calculate="calculateNewMeasurement"
          :measurementTimestamps="props.plate.details.measurement_timestamps" />
      </q-card-section>
    </q-card>
  </q-dialog>
</template>

<style scoped lang="sass">

.inner
  border-radius: 50%
  line-height: 18px
  margin: auto
  text-align: center
  width: 22px
  height: 22px
  min-width: 22px
  min-height: 22px
  max-width: 22px
  max-height: 22px


table
  border-spacing: 5px
  border: 1px solid #bbb
  border-radius: 12px
  padding: 4px 8px 8px 4px
  overflow: hidden


td

  width: 24px
  height: 24px
  min-width: 24px
  min-height: 24px
  max-width: 24px
  max-height: 24px
  text-align: center
  vertical-align: middle
  font-size: 8px
td
  a
    text-decoration: none
td.circle
  border: 1px solid #4c4c4c
  border-radius: 12px



td:hover
  background-color: #bbb
  cursor: pointer
ul
  list-style: none
  padding-left: 0
.tooltip
  table, tr, td, td:hover
    font-size: 12px
    border: unset
    border-radius: unset
    border-spacing: 2px
    padding: 0px
    margin: 0px
    overflow: unset
    text-align: left
    vertical-align: top
  td
    white-space: nowrap
    width: unset
    height: unset
    min-width: unset
    min-height: unset
    max-width: unset
    max-height: unset
  & > b
    font-size: 12px
    margin-top: 10px

select
  max-width: 300px
.measurement
  border-radius: 5px


.calculator_dialog
  min-width: 800px

table.smaller
  border-spacing: 0px
  border-collapse: collapse
  border: 0px solid #eee
  border-radius: unset
  padding: unset
  overflow: hidden

.smaller td
  width: 15px
  height: 15px
  min-width: 15px
  min-height: 15px
  max-width: 15px
  max-height: 15px
  border-width: 1px
  border-color: #eee
  border-radius: 0
  text-align: center
  vertical-align: middle
  font-size: 8px
  padding: 0

.smaller th
  font-size: 8px
</style>
