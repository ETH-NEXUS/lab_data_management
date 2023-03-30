<script setup lang="ts">
import {computed, defineProps, defineEmits, PropType, ref, onMounted} from 'vue'
import {Plate, Well, LegendColor, MeasurementInfo, MinMax, Measurement} from './models'
import {positionFromRowCol} from '../helpers/plate'
import {palettes, percentageToHsl} from 'components/helpers'
import {storeToRefs} from 'pinia'
import {useSettingsStore} from 'stores/settings'
import {useI18n} from 'vue-i18n'

const {t} = useI18n()

const {platePage} = storeToRefs(useSettingsStore())

const props = defineProps({
  plate: {
    type: Object as PropType<Plate>,
    required: true,
  },
})

//computed message saying that the combination of label and timestamp was not found

const notFoundMessage = computed(() => {
  const measurement = props.plate.wells![0].measurements?.find(
    (m: Measurement) =>
      m.feature.abbrev === selectedLabel.value && m.measurement_timestamp === selectedTimestamp.value
  )

  if (!measurement) {
    return `No data found for ${selectedLabel.value} at ${selectedTimestamp.value}`
  }
  return ''
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

const selectedMeasurement = ref<MeasurementInfo | undefined>()
const labels = ref<Array<string>>([])
const timestamps = ref<Array<string>>([])
const measurementOptions = ref<Array<MeasurementInfo>>([])
const selectedLabel = ref<string>('')
const selectedTimestamp = ref<string>('')

const legendColors = computed(() => {
  return createColorLegend()
})

const emit = defineEmits(['well-selected'])

const alpha = Array.from(Array(26))
  .map((e, i) => i + 65)
  .map(x => String.fromCharCode(x))

const byPosition = (position: number) => {
  return props.plate.wells?.find(w => w.position === position)
}

const wells = computed(() => {
  const _wells: Array<Array<Well | undefined>> = Array.from(Array(props.plate.dimension.rows), () =>
    new Array(props.plate.dimension.cols).fill(undefined)
  )
  for (const row of Array(props.plate.dimension.rows).keys()) {
    for (const col of Array(props.plate.dimension.cols).keys()) {
      _wells[row][col] = byPosition(row * props.plate.dimension.cols + col)
      // _wells[row][col]?.measurements.push({
      //   value: (row * props.plate.dimension.cols + col) / 352,
      //   name: '',
      //   abbrev: '',
      //   unit: '',
      // })
    }
  }
  return _wells
})

const percentage = (well: Well | undefined) => {
  if (selectedLabel.value && selectedTimestamp.value && well?.measurements) {
    const min = props.plate.min_max.find(
      (item: MinMax) => item.feature === selectedLabel.value && item.timestamp === selectedTimestamp.value
    )?.min_all_types

    const max = props.plate.min_max.find(
      (item: MinMax) => item.feature === selectedLabel.value && item.timestamp === selectedTimestamp.value
    )?.max_all_types

    const value = well?.measurements.filter(
      m => m.feature.abbrev === selectedLabel.value && m.measurement_timestamp === selectedTimestamp.value
    )[0]?.value

    if (value && min && max) {
      return (value - min) / (max - min)
    }
  }
  return 0
}

const heatmapColor = (well: Well | undefined) => {
  const {from, to} = platePage.value.heatmapPalette.value
  return percentageToHsl(percentage(well), from, to)
}

const typeColor_map: {[key: string]: string} = {
  P: 'rgb(198, 223, 168)',
  N: 'rgb(253, 204, 134)',
  C: 'rgb(177, 190, 197)',
}

const typeColor = (well: Well | undefined) => {
  if (well) {
    const type = well.type.substring(0, 1).toLocaleUpperCase()
    if (type in typeColor_map) {
      return typeColor_map[type]
    }
  }
  return 'transparent'
}

const createColorLegend = () => {
  const legend: LegendColor[] = []
  if (selectedLabel.value && selectedTimestamp.value) {
    const min = props.plate.min_max.find(
      (item: MinMax) => item.feature === selectedLabel.value && item.timestamp === selectedTimestamp.value
    )?.min_all_types

    const max = props.plate.min_max.find(
      (item: MinMax) => item.feature === selectedLabel.value && item.timestamp === selectedTimestamp.value
    )?.max_all_types

    const numberOfSteps = 10
    if (min && max) {
      const step: number = (max - min) / numberOfSteps
      for (let i = 0; i <= numberOfSteps; i++) {
        const value: number = min + i * step
        const color = percentageToHsl(
          (value - min) / (max - min),
          platePage.value.heatmapPalette.value.from,
          platePage.value.heatmapPalette.value.to
        )
        legend.push({value, color})
      }
    }

    return legend.reverse()
  }
}

onMounted(() => {
  if (props.plate.measurements && props.plate.measurements.length > 0) {
    selectedMeasurement.value = props.plate.measurements[0]

    const _labels = props.plate?.measurements.map(m => m.feature)
    labels.value = [...new Set(_labels)]
    selectedLabel.value = props.plate.measurements[0].feature

    const _timestamps = props.plate?.measurements.map(m => m.measurement_timestamp)
    timestamps.value = [...new Set(_timestamps)]
    selectedTimestamp.value = props.plate.measurements[0].measurement_timestamp

    measurementOptions.value = props.plate.measurements.map(m => {
      return {
        feature: m.feature,
        label: m.feature,
        value: m.feature,
        measurement_timestamp: m.measurement_timestamp,
      }
    })
  }
})
</script>

<template>
  <div class="row">
    <table v-if="props.plate">
      <tr>
        <th />
        <th :key="`cols${col}`" v-for="(_, col) of props.plate.dimension.cols">
          {{ col + 1 }}
        </th>
      </tr>
      <tr :key="`row${row}`" v-for="(_, row) of props.plate.dimension.rows">
        <th>
          {{ alpha[row] }}
        </th>

        <td
          :style="{
            backgroundColor:
              platePage.showHeatmap && selectedLabel
                ? heatmapColor(wells[row][col])
                : plate.template || (platePage.wellContent === 'type' && !platePage.showHeatmap)
                ? typeColor(wells[row][col])
                : 'transparent',
          }"
          :key="`cols${col}`"
          v-for="(_, col) of props.plate.dimension.cols"
          @click="
            emit('well-selected', {
              well: wells[row][col],
              position: positionFromRowCol(row, col, props.plate.dimension),
            })
          ">
          <a v-if="wells[row][col]" :class="{'bg-warning': wells[row][col]!.status}">
            {{ wells[row][col]![platePage.wellContent] }}
          </a>
          <q-tooltip
            class="tooltip"
            v-if="wells[row][col]"
            anchor="top middle"
            self="bottom middle"
            :offset="[5, 5]">
            <b>{{ wells[row][col]?.hr_position }}</b>
            ({{ wells[row][col]?.type }})
            <small>({{ wells[row][col]?.position }})</small>
            <hr />
            <b v-if="wells[row][col]?.compounds" class="q-mt-md">{{ t('label.compounds') }}</b>
            <table v-if="wells[row][col]?.compounds">
              <tr v-for="compound in wells[row][col]?.compounds" :key="compound.identifier">
                <td>
                  <b>{{ compound.name }}</b>
                </td>
                <td>{{ compound.identifier }}</td>
              </tr>
            </table>
            <b v-if="wells[row][col]?.measurements" class="q-mt-md">{{ t('label.measurements') }}</b>
            <table v-if="wells[row][col]?.measurements">
              <tr v-for="measurement in wells[row][col]?.measurements" :key="measurement.feature.name">
                <td>
                  <b>{{ measurement.feature.abbrev }}</b>
                </td>
                <td>{{ measurement.value }}</td>
                <td>{{ measurement.feature.unit }}</td>
              </tr>
            </table>
          </q-tooltip>
        </td>
      </tr>
    </table>
    <div v-if="platePage.showHeatmap && labels.length > 0 && legendColors" class="q-my-md q-ml-md">
      <div
        class="legendItem"
        v-for="(color, idx) in legendColors"
        :key="color.value + idx"
        :style="{backgroundColor: color.color}">
        <span class="legendLabel">{{ [0, 5, 10].includes(idx) ? color.value.toFixed(0) : ' ' }}</span>
      </div>
    </div>
  </div>
  <div v-if="!plate.template && props.plate.z_primes" class="row q-mt-sm">
    <div class="col-12 q-mr-md">
      <b>{{ t('label.z_prime') }}:</b>
      <span
        v-for="item in props.plate.z_primes"
        :key="item.feature"
        :class="{
          'z-prime': true,
          'q-mx-xs': true,
          'q-pa-xs': true,
          'bg-green-3': item.z_prime >= 0.5 && item.z_prime <= 1,
          'bg-orange-3': item.z_prime >= 0 && item.z_prime < 0.5,
          'bg-red-3': item.z_prime < 0,
        }">
        {{ item.z_prime }} ({{ item.feature }}, {{ item.timestamp }})
      </span>
    </div>
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
    <div v-if="measurementOptions.length > 0" class="col-4">
      <q-checkbox v-model="platePage.showHeatmap" :label="t('label.show_heatmap')"></q-checkbox>
      <div class="col-12" v-if="platePage.showHeatmap && measurementOptions.length">
        <q-select
          v-model="selectedLabel"
          :options="labels"
          :label="t('message.select_label')"
          v-if="labels.length > 1"></q-select>
        <q-select
          v-if="timestamps.length > 1"
          v-model="selectedTimestamp"
          :options="timestamps"
          :label="t('message.select_timestamp')"></q-select>
      </div>
      <q-banner inline-actions class="text-white bg-red" v-if="notFoundMessage && platePage.showHeatmap">
        {{ notFoundMessage }}
      </q-banner>

      <div class="col-12" v-if="platePage.showHeatmap && measurementOptions.length > 0">
        <q-select
          v-model="platePage.heatmapPalette"
          :options="palettes"
          :label="t('label.select_color_palette')"></q-select>
      </div>
    </div>
  </div>
</template>

<style scoped lang="sass">
table
  border-spacing: 5px
  border: 1px solid #bbb
  border-radius: 12px
  padding: 4px 8px 8px 4px
  overflow: hidden
td
  border: 1px solid #4c4c4c
  border-radius: 12px
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
.z-prime
  border-radius: 5px

.legendItem
  position: relative
  width: 50px
  height: 20px

.legendLabel
  position: absolute
  left: 53px
</style>

<!-- :style="
            'background-color: ' + percentageToHsl(wells[row][col]?.measurements?[0]?.value || -1, 120, 0)
          " -->
