<script setup lang="ts">
import {computed, defineProps, defineEmits, PropType, ref, onMounted} from 'vue'
import {Plate, WellDetails, LegendColor} from './models'
import {positionFromRowCol} from '../helpers/plate'
import {palettes, percentageToHsl} from 'components/helpers'
import {storeToRefs} from 'pinia'
import {useSettingsStore} from 'stores/settings'
import {useI18n} from 'vue-i18n'

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

const $q = useQuasar()

const selectedMeasurement = ref<string | null>(null)
const measurementOptions = ref<string[] | null>(null)
const selectedTimestampIdx = ref<number>(0)

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
</script>

<template>
  <div class="row no-wrap">
    <table v-if="props.plate">
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
            backgroundColor: heatmapColor(wells[row][col]),
          }"
          :key="`cols${col}`"
          v-for="(_, col) of props.plate.dimension.cols">
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

    <div class="q-my-md q-ml-md">
      <div
        class="legendItem"
        v-for="(color, idx) in legendColors"
        :key="color.value + idx"
        :style="{backgroundColor: color.color}">
        <span class="legendLabel">{{ [0, 5, 10].includes(idx) ? color.value.toFixed(3) : ' ' }}</span>
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
