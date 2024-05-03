<script setup lang="ts">
import {computed, defineEmits, defineProps, PropType} from 'vue'
import {positionFromRowCol} from 'src/helpers/plate'
import WellTooltip from 'components/wells/WellTooltip.vue'
import {storeToRefs} from 'pinia'
import {useSettingsStore} from 'stores/settings'
import {Plate, WellDetails} from 'components/models'
import {percentageToHsl} from 'components/helpers'

const {platePage} = storeToRefs(useSettingsStore())

const props = defineProps({
  plate: {
    type: Object as PropType<Plate>,
    required: true,
  },
  plateIndex: {
    type: Number,
    default: 0,
  },
  min: {
    type: Number,
    required: true,
  },
  max: {
    type: Number,
    required: true,
  },
  selectedMeasurement: {
    type: String as () => string | null,
    default: null,
  },
  selectedTimestampIdx: {
    type: Number,
    default: null,
  },
})
const emit = defineEmits(['well-selected'])

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
const heatmapColor = (well: WellDetails | undefined) => {
  const {from, to} = platePage.value.heatmapPalette.value

  return percentageToHsl(percentage(well), from, to)
}

const typeColor_map: {[key: string]: string} = {
  P: 'rgb(107,142,35)',
  N: 'rgb(255, 99, 71)',
  C: 'rgb(255,255,255)',
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

const measurement = (well: WellDetails) => {
  if (props.selectedMeasurement && props.selectedTimestampIdx !== null && well.measurements) {
    // it should be !== null, because the index can be 0
    if (props.selectedMeasurement in well.measurements) {
      if (well.measurements[props.selectedMeasurement].length > props.selectedTimestampIdx) {
        return well.measurements[props.selectedMeasurement][props.selectedTimestampIdx]
      }
    }
  }
  return null
}

const percentage = (well: WellDetails | undefined) => {
  console.log(well)
  if (well) {
    const value = measurement(well)
    if (value != null) {
      return (value - props.min) / (props.max - props.min)
    }
  }
  return -1
}
</script>

<template>
  <table :class="platePage.smallerMapView ? 'smaller' : ''">
    <tr>
      <th />
      <th :key="`cols${col}-${plateIndex}`" v-for="(_, col) of props.plate.dimension.cols">
        {{ col + 1 }}
      </th>
    </tr>
    <tr :key="`row${row}-${plateIndex}`" v-for="(_, row) of props.plate.dimension.rows">
      <th>
        {{ posToAlphaChar(row + 1) }}
      </th>
      <td
        :class="platePage.squareCompoundType && !platePage.smallerMapView ? 'square' : 'circle'"
        :style="{
          backgroundColor:
            platePage.squareCompoundType && !platePage.smallerMapView
              ? typeColor(wells[row][col])
              : platePage.showHeatmap && selectedMeasurement
              ? heatmapColor(wells[row][col])
              : plate.template || (platePage.wellContent === 'type' && !platePage.showHeatmap)
              ? typeColor(wells[row][col])
              : 'transparent',
        }"
        :key="`cols_${col}-${plateIndex}`"
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
            :col="Number(col)"
            :row="Number(row)" />
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
            :col="Number(col)"
            :row="Number(row)" />
        </a>
      </td>
    </tr>
  </table>
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
  border: 1px solid #4c4c4c


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
