<script setup lang="ts">
import { computed } from 'vue'
import type { Plate, WellDetails, WellInfo } from '~/types/lab'
import { usePlateViewStore } from '~/stores/plateView'
import WellTooltip from '~/components/wells/WellTooltip.vue'
import { percentageToHsl } from '~/utils/plateHeatmap'

type Props = {
  plate: Plate
  plateIndex?: number
  min: number
  max: number
}

const props = withDefaults(defineProps<Props>(), {
  plateIndex: 0,
})

const emit = defineEmits<{
  (e: 'well-selected', payload: WellInfo): void
}>()

const platePage = usePlateViewStore()

const byPosition = (position: number) => {
  return props.plate.wells?.find((well) => well.position === position)
}

const wells = computed(() => {
  const tableWells: Array<Array<WellDetails | undefined>> = Array.from(Array(props.plate.dimension.rows), () =>
    new Array(props.plate.dimension.cols).fill(undefined),
  )
  for (const row of Array(props.plate.dimension.rows).keys()) {
    for (const col of Array(props.plate.dimension.cols).keys()) {
      tableWells[row]![col] = byPosition(row * props.plate.dimension.cols + col)
    }
  }

  return tableWells
})

const asciiUppercase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

const posToAlphaChar = (pos: number) => {
  let letter = ''
  while (pos > 0) {
    const remainder = (pos - 1) % 26
    pos = Math.floor((pos - 1) / 26)
    letter = asciiUppercase[remainder] + letter
  }
  return letter
}

const heatmapColor = (well: WellDetails | undefined) => {
  const { from, to } = platePage.heatmapPalette.value
  return percentageToHsl(percentage(well), from, to, platePage.heatmapPalette.label)
}

const typeColorMap: { [key: string]: string } = {
  P: 'rgb(107,142,35)',
  N: 'rgb(255, 99, 71)',
  C: 'rgb(255,255,255)',
}

const typeColor = (well: WellDetails | undefined) => {
  if (well) {
    const type = well.type.substring(0, 1).toLocaleUpperCase()
    if (type in typeColorMap) {
      return typeColorMap[type] as string
    }
  }
  return 'transparent'
}

const measurement = (well: WellDetails) => {
  const selectedMeasurement = platePage.selectedMeasurement
  const selectedTimestampIdx = platePage.selectedTimestampIdx

  if (!selectedMeasurement || !well.measurements) {
    return null
  }

  // it should be !== null, because the index can be 0
  if (!(selectedMeasurement in well.measurements)) {
    return null
  }

  const series = well.measurements[selectedMeasurement]
  if (!series || series.length <= selectedTimestampIdx) {
    return null
  }

  return series[selectedTimestampIdx] ?? null
}

const percentage = (well: WellDetails | undefined) => {
  if (well) {
    const value = measurement(well)
    if (value != null) {
      return (value - props.min) / (props.max - props.min)
    }
  }
  return -1
}

const positionFromRowCol = (row: number, col: number) => row * props.plate.dimension.cols + col
</script>

<template>
  <table :class="platePage.smallerMapView ? 'smaller' : ''">
    <tr>
      <th />
      <th v-for="(_colItem, col) of props.plate.dimension.cols" :key="`cols${col}-${plateIndex}`">
        {{ col + 1 }}
      </th>
    </tr>
    <tr v-for="(_rowItem, row) of props.plate.dimension.rows" :key="`row${row}-${plateIndex}`">
      <th>
        {{ posToAlphaChar(row + 1) }}
      </th>
      <td
        v-for="(_, col) of props.plate.dimension.cols"
        :key="`cols_${col}-${plateIndex}`"
        :class="platePage.squareCompoundType && !platePage.smallerMapView ? 'square' : 'circle'"
        :style="{
          backgroundColor:
            platePage.squareCompoundType && !platePage.smallerMapView
              ? typeColor(wells[row]?.[col])
              : platePage.showHeatmap && platePage.selectedMeasurement
                ? heatmapColor(wells[row]?.[col])
                : props.plate.template || (platePage.wellContent === 'type' && !platePage.showHeatmap)
                  ? typeColor(wells[row]?.[col])
                  : 'transparent',
        }"
        @click="
          emit('well-selected', {
            well: wells[row]?.[col] as WellDetails,
            position: positionFromRowCol(row, col),
          })
        "
      >
        <div
          v-if="platePage.squareCompoundType && !platePage.smallerMapView"
          :class="platePage.smallerMapView ? 'innerSmaller' : 'inner'"
          :style="{
            backgroundColor:
              platePage.showHeatmap && platePage.selectedMeasurement
                ? heatmapColor(wells[row]?.[col])
                : props.plate.template || (platePage.wellContent === 'type' && !platePage.showHeatmap)
                  ? typeColor(wells[row]?.[col])
                  : 'transparent',
          }"
        >
          <WellTooltip v-if="wells[row]?.[col]" :well="wells[row]?.[col]" :col="Number(col)" :row="Number(row)">
            <a :class="{ 'bg-warning': wells[row]?.[col]?.status }">
              {{ platePage.smallerMapView ? '&nbsp;&nbsp;&nbsp;' : wells[row]?.[col]?.[platePage.wellContent] }}
            </a>
          </WellTooltip>
        </div>
        <WellTooltip
          v-if="wells[row]?.[col] && !platePage.squareCompoundType"
          :well="wells[row]?.[col]"
          :col="Number(col)"
          :row="Number(row)"
        >
          <a :class="{ 'bg-warning': wells[row]?.[col]?.status }">
            {{ platePage.smallerMapView ? '&nbsp;&nbsp;&nbsp;' : wells[row]?.[col]?.[platePage.wellContent] }}
          </a>
        </WellTooltip>
      </td>
    </tr>
  </table>
</template>

<style scoped>
.inner {
  border-radius: 50%;
  line-height: 18px;
  margin: auto;
  text-align: center;
  width: 22px;
  height: 22px;
  min-width: 22px;
  min-height: 22px;
  max-width: 22px;
  max-height: 22px;
  border: 1px solid #4c4c4c;
}

table {
  border-collapse: separate;
  border-spacing: 4px;
  background-color: #fff;
  border: 1px solid #bbb;
  border-radius: 12px;
  padding: 4px 8px 8px 4px;
  overflow: hidden;
}

th {
  font-size: 14px;
}

td {
  width: 24px;
  height: 24px;
  min-width: 24px;
  min-height: 24px;
  max-width: 24px;
  max-height: 24px;
  text-align: center;
  vertical-align: middle;
  font-size: 8px;
}

td a {
  text-decoration: none;
}

.bg-warning {
  background-color: #ffc107;
  color: #000;
}

td.circle {
  border: 1px solid #4c4c4c;
  border-radius: 12px;
}

td:hover {
  background-color: #bbb;
  cursor: pointer;
}

ul {
  list-style: none;
  padding-left: 0;
}

.tooltip table,
.tooltip tr,
.tooltip td,
.tooltip td:hover {
  font-size: 12px;
  border: unset;
  border-radius: unset;
  border-spacing: 2px;
  padding: 0;
  margin: 0;
  overflow: unset;
  text-align: left;
  vertical-align: top;
}

.tooltip td {
  white-space: nowrap;
  width: unset;
  height: unset;
  min-width: unset;
  min-height: unset;
  max-width: unset;
  max-height: unset;
}

.tooltip > b {
  font-size: 12px;
  margin-top: 10px;
}

select {
  max-width: 300px;
}

.measurement {
  border-radius: 5px;
}

.calculator_dialog {
  min-width: 800px;
}

table.smaller {
  border-spacing: 0;
  border-collapse: collapse;
  border: 0 solid #eee;
  border-radius: unset;
  padding: unset;
  overflow: hidden;
}

.smaller td {
  width: 15px;
  height: 15px;
  min-width: 15px;
  min-height: 15px;
  max-width: 15px;
  max-height: 15px;
  border-width: 1px;
  border-color: #eee;
  border-radius: 0;
  text-align: center;
  vertical-align: middle;
  font-size: 8px;
  padding: 0;
}

.smaller th {
  font-size: 9px;
}
</style>
