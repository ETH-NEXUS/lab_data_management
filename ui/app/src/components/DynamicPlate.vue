<script setup lang="ts">
import {computed, defineProps, defineEmits, PropType, ref, onMounted} from 'vue'
import {Plate, WellDetails, LegendColor} from './models'
import {positionFromRowCol} from '../helpers/plate'
import {palettes, percentageToHsl} from 'components/helpers'
import {storeToRefs} from 'pinia'
import {useSettingsStore} from 'stores/settings'
import {useI18n} from 'vue-i18n'
import MeasurementCalculator from 'components/MeasurementCalculator.vue'
import {useProjectStore} from 'stores/project'
import bus from 'src/eventBus'
import {useQuasar} from 'quasar'
import WellTooltip from 'components/WellTooltip.vue'

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
  // const combined_labels: string[] = []

  // for (let i = 0; i < props.plate.details.measurement_labels.length; ++i) {
  //   if (props.plate.details.measurement_timestamps && props.plate.details.measurement_timestamps.length > 1) {
  //     for (let j = 0; j < props.plate.details.measurement_timestamps.length; ++j) {
  //       combined_labels.push(
  //         `${props.plate.details.measurement_labels[i]} (${props.plate.details.measurement_timestamps[j]})`
  //       )
  //     }
  //   } else {
  //     combined_labels.push(props.plate.details.measurement_labels[i])
  //   }
  // }

  // return combined_labels
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

const legendColors = computed(() => {
  return createColorLegend()
})

const emit = defineEmits(['well-selected', 'refresh'])

const byPosition = (position: number) => {
  return props.plate.wells?.find(w => w.position === position)
}

const wells = computed(() => {
  const _wells: Array<Array<WellDetails | undefined>> = Array.from(Array(props.plate.dimension.rows), () =>
    new Array(props.plate.dimension.cols).fill(undefined)
  )
  for (const row of Array(props.plate.dimension.rows).keys()) {
    for (const col of Array(props.plate.dimension.cols).keys()) {
      _wells[row][col] = byPosition(row * props.plate.dimension.cols + col)
    }
  }
  return _wells
})

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

const ascii_uppercase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

const posToAlphaChar = (pos: number) => {
  let letter = ''
  while (pos > 0) {
    let remainder = (pos - 1) % 26
    pos = Math.floor((pos - 1) / 26)
    letter = ascii_uppercase[remainder] + letter
  }
  return letter
}

const min = computed<number>(() => {
  if (!selectedMeasurement.value) return 0
  return props.plate.details.overall_stats[selectedMeasurement.value].min[selectedTimestampIdx.value]
})

const max = computed<number>(() => {
  if (!selectedMeasurement.value) return 0
  return props.plate.details.overall_stats[selectedMeasurement.value].max[selectedTimestampIdx.value]
})

const measurement = (well: WellDetails) => {
  if (selectedMeasurement.value) {
    if (selectedMeasurement.value in well.measurements) {
      if (well.measurements[selectedMeasurement.value].length > selectedTimestampIdx.value) {
        return well.measurements[selectedMeasurement.value][selectedTimestampIdx.value]
      }
    }
  }
  return null
}

const percentage = (well: WellDetails | undefined) => {
  if (well) {
    const value = measurement(well)
    if (value) {
      return (value - min.value) / (max.value - min.value)
    }
  }
  return 0
}

const heatmapColor = (well: WellDetails | undefined) => {
  const {from, to} = platePage.value.heatmapPalette.value
  return percentageToHsl(percentage(well), from, to)
}

const typeColor_map: {[key: string]: string} = {
  P: 'rgb(107,142,35)',
  N: 'rgb(255, 99, 71)',
  C: 'rgb(0,191,255)',
}

const typeColor = (well: WellDetails | undefined) => {
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
  if (selectedMeasurement.value) {
    const numberOfSteps = 10

    const step: number = (max.value - min.value) / numberOfSteps
    for (let i = 0; i <= numberOfSteps; i++) {
      const value: number = min.value + i * step
      const color = percentageToHsl(
        (value - min.value) / (max.value - min.value),
        platePage.value.heatmapPalette.value.from,
        platePage.value.heatmapPalette.value.to
      )
      legend.push({value, color})
    }

    return legend.reverse()
  }
}

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
  <q-toggle
    class="q-mb-md"
    size="sm"
    checked-icon="check"
    v-model="platePage.smallerMapView"
    :label="t('label.smaller_map_view')"
    right-label
    color="secondary"></q-toggle>
  <q-toggle
    class="q-mb-md"
    size="sm"
    checked-icon="check"
    v-model="platePage.squareCompoundType"
    :label="t('label.show_type_as_square')"
    right-label
    color="secondary"></q-toggle>
  <div class="row no-wrap">
    <table v-if="props.plate" :class="platePage.smallerMapView ? 'smaller' : ''">
      <tr>
        <th />
        <th :key="`cols${col}`" v-for="(_, col) of props.plate.dimension.cols">
          {{ col + 1 }}
        </th>
      </tr>
      <tr :key="`row${row}`" v-for="(_, row) of props.plate.dimension.rows">
        <th>
          {{ posToAlphaChar(row + 1) }}
        </th>

        <td
          :class="platePage.squareCompoundType ? 'square' : 'circle'"
          :style="{
            backgroundColor: platePage.squareCompoundType
              ? typeColor(wells[row][col])
              : platePage.showHeatmap && selectedMeasurement
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
          <div
            v-if="platePage.squareCompoundType && !platePage.smallerMapView"
            :class="platePage.smallerMapView ? 'innerSmaller' : 'inner'"
            :style="{
              backgroundColor:
                platePage.showHeatmap && selectedMeasurement
                  ? heatmapColor(wells[row][col])
                  : plate.template || (platePage.wellContent === 'type' && !platePage.showHeatmap)
                  ? typeColor(wells[row][col])
                  : 'transparent',
            }">
            <a v-if="wells[row][col]" :class="{'bg-warning': wells[row][col]!.status}">
              {{ platePage.smallerMapView ? '&nbsp;&nbsp;&nbsp;' : wells[row][col]![platePage.wellContent] }}
            </a>
            <WellTooltip
              v-if="wells[row][col]"
              :well="wells[row][col]"
              :selected-timestamp-idx="selectedTimestampIdx"
              :selected-measurement="selectedMeasurement"
              :col="col"
              :row="row" />
          </div>
          <a
            v-if="wells[row][col] && !platePage.squareCompoundType"
            :class="{'bg-warning': wells[row][col]!.status}">
            {{ platePage.smallerMapView ? '&nbsp;&nbsp;&nbsp;' : wells[row][col]![platePage.wellContent] }}
            <WellTooltip
              v-if="wells[row][col] && !platePage.squareCompoundType"
              :well="wells[row][col]"
              :selected-timestamp-idx="selectedTimestampIdx"
              :selected-measurement="selectedMeasurement"
              :col="col"
              :row="row" />
          </a>
        </td>
      </tr>
    </table>

    <div v-if="platePage.showHeatmap && selectedMeasurement && legendColors" class="q-my-md q-ml-md">
      <div
        class="legendItem"
        v-for="(color, idx) in legendColors"
        :key="color.value + idx"
        :style="{backgroundColor: color.color}">
        <span class="legendLabel">{{ [0, 5, 10].includes(idx) ? color.value.toFixed(3) : ' ' }}</span>
      </div>
    </div>
  </div>
  <div class="row">
    <div v-if="z_prime" class="q-mt-sm">
      <b>{{ t('label.z_prime') }}:</b>
      <span
        :class="{
          'measurement': true,
          'q-mx-xs': true,
          'q-pa-xs': true,
          'bg-green-3': z_prime >= 0.5 && z_prime <= 1,
          'bg-orange-3': z_prime >= 0 && z_prime < 0.5,
          'bg-red-3': z_prime < 0,
        }">
        {{ z_prime.toFixed(2) }}
      </span>
    </div>
    <div v-if="ssmd" class="q-mt-sm">
      <b>{{ t('label.ssmd') }}:</b>
      <span
        :class="{
          'measurement': true,
          'q-mx-xs': true,
          'q-pa-xs': true,
          'bg-green-3': ssmd >= 12,
          'bg-orange-3': ssmd >= 6 && ssmd < 12,
          'bg-red-3': ssmd < 6,
        }">
        {{ ssmd.toFixed(2) }}
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

    <div v-if="measurementOptions" class="col-4">
      <q-checkbox
        v-model="platePage.showHeatmap"
        :label="t('label.show_heatmap')"
        class="q-mr-md"></q-checkbox>

      <div class="col-12" v-if="platePage.showHeatmap && measurementOptions.length">
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
        <MeasurementCalculator :labels="combinedLabels" @calculate="calculateNewMeasurement" />
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

.innerSmaller
  line-height: 15px
  width: 15px
  height: 15px
  min-width: 15px
  min-height: 15px
  max-width: 15px
  max-height: 15px

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

.legendItem
  position: relative
  width: 30px
  height: 15px

.legendLabel
  position: absolute
  left: 33px
  font-size: 10px


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

<!-- :style="
            'background-color: ' + percentageToHsl(wells[row][col]?.measurements?[0]?.value || -1, 120, 0)
          " -->
