<script setup lang="ts">
import {computed, defineProps, defineEmits, PropType, ref, onMounted} from 'vue'
import {Plate, WellDetails, LegendColor, SelectOption} from './models'
import {positionFromRowCol} from '../helpers/plate'
import {palettes, percentageToHsl} from 'components/helpers'
import {storeToRefs} from 'pinia'
import {useSettingsStore} from 'stores/settings'
import {useI18n} from 'vue-i18n'
import MeasurementCalculator from 'components/MeasurementCalculator.vue'
import {useProjectStore} from 'stores/project'

const {t} = useI18n()

const {platePage} = storeToRefs(useSettingsStore())
const projectStore = useProjectStore()

const props = defineProps({
  plate: {
    type: Object as PropType<Plate>,
    required: true,
  },
})

//computed message saying that the combination of label and timestamp was not found

// const notFoundMessage = computed(() => {
//   if (!props.plate.wells || !(selectedMeasurement.value in props.plate.wells[0].measurements)) {
//     return `No data found for ${selectedMeasurement.value} at ${selectedTimestampIdx.value}`
//   }
//   return ''
// })

const combinedLabels = computed<string[]>(() => {
  const combined_labels: string[] = []

  for (let i = 0; i < props.plate.details.measurement_labels.length; ++i) {
    if (props.plate.details.measurement_timestamps && props.plate.details.measurement_timestamps.length > 1) {
      for (let j = 0; j < props.plate.details.measurement_timestamps.length; ++j) {
        combined_labels.push(
          `${props.plate.details.measurement_labels[i]} (${props.plate.details.measurement_timestamps[j]})`
        )
      }
    } else {
      combined_labels.push(props.plate.details.measurement_labels[i])
    }
  }

  return combined_labels
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

const selectedMeasurement = ref<string | null>(null)
const measurementOptions = ref<string[]>([])

const selectedTimestampIdx = ref<number>(0)
const timestampOptions = ref<Array<SelectOption<number>>>([])

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
  P: 'rgb(198, 223, 168)',
  N: 'rgb(253, 204, 134)',
  C: 'rgb(177, 190, 197)',
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
  if (props.plate.details.measurement_timestamps && props.plate.details.measurement_timestamps.length > 0) {
    timestampOptions.value = props.plate.details.measurement_timestamps.map((ts, idx) => {
      return {label: ts, value: idx}
    })
    selectedTimestampIdx.value = 0
  }
})

const calculateNewMeasurement = async (expression: string, newLabel: string) => {
  await projectStore.addNewMeasurement(props.plate.id, expression, newLabel)
  openCalculator.value = false

  emit('refresh')
}
</script>

<template>
  <div class="row">
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
          :style="{
            backgroundColor:
              platePage.showHeatmap && selectedMeasurement
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
            {{ platePage.smallerMapView ? '&nbsp;&nbsp;&nbsp;' : wells[row][col]![platePage.wellContent] }}
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
              <tr v-for="compound in wells[row][col]?.compounds" :key="compound">
                <b>{{ compound }}</b>
              </tr>
            </table>
            <b v-if="selectedMeasurement" class="q-mt-md">{{ t('label.measurements') }}</b>
            <br />
            <span v-if="selectedMeasurement" class="q-mt-md">
              {{ selectedMeasurement }}: {{ measurement(wells[row][col]!) }}
            </span>
          </q-tooltip>
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
  <div v-if="z_prime" class="row q-mt-sm">
    <div class="col-12 q-mr-md">
      <b>{{ t('label.z_prime') }}:</b>
      <span
        :class="{
          'z-prime': true,
          'q-mx-xs': true,
          'q-pa-xs': true,
          'bg-green-3': z_prime >= 0.5 && z_prime <= 1,
          'bg-orange-3': z_prime >= 0 && z_prime < 0.5,
          'bg-red-3': z_prime < 0,
        }">
        {{ z_prime }} ({{ selectedMeasurement }},
        {{ timestampOptions.find(tso => tso.value === selectedTimestampIdx)?.label }})
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
      <q-checkbox
        v-model="platePage.smallerMapView"
        :label="t('label.smaller_map_view')"
        class="q-mt-md"></q-checkbox>
    </div>

    <div v-if="measurementOptions.length > 0" class="col-4">
      <q-checkbox
        v-model="platePage.showHeatmap"
        :label="t('label.show_heatmap')"
        class="q-mr-md"></q-checkbox>

      <q-btn
        v-if="platePage.showHeatmap"
        class="q-my-md"
        :label="t('action.calculate_measurement')"
        icon="calculate"
        color="primary"
        @click="openCalculator = true" />

      <div class="col-12" v-if="platePage.showHeatmap && measurementOptions.length">
        <q-select
          v-model="selectedMeasurement"
          :options="measurementOptions"
          :label="t('message.select_label')"
          v-if="measurementOptions.length > 1"></q-select>
        <q-select
          v-if="timestampOptions.length > 1"
          v-model="selectedTimestampIdx"
          :options="timestampOptions"
          :label="t('message.select_timestamp')"
          emit-value
          map-options></q-select>
      </div>
      <!-- <q-banner inline-actions class="text-white bg-red" v-if="notFoundMessage && platePage.showHeatmap">
        {{ notFoundMessage }}
      </q-banner> -->

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
