<script setup lang="ts">
import {computed, defineProps, defineEmits, PropType} from 'vue'
import {Plate, Well} from './models'
import {positionFromRowCol} from '../helpers/plate'

const props = defineProps({
  plate: {
    type: Object as PropType<Plate>,
    required: true,
  },
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

const percentageToHsl = (percentage: number, hue0: number, hue1: number) => {
  // if percentage is not given (-1) we return a transparent color
  if (percentage === -1) {
    return 'rgba(255,255,255,0)'
  }
  var hue = percentage * (hue1 - hue0) + hue0
  return 'hsl(' + hue + ', 100%, 50%)'
}
</script>

<template>
  <div>
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
          :key="`cols${col}`"
          v-for="(_, col) of props.plate.dimension.cols"
          @click="
            emit('well-selected', {
              well: wells[row][col],
              position: positionFromRowCol(row, col, props.plate.dimension),
            })
          ">
          <a v-if="wells[row][col]">
            {{ wells[row][col]!.hr_position }}
          </a>
          <q-tooltip
            v-if="wells[row][col] && wells[row][col]?.compounds"
            anchor="top middle"
            self="bottom middle"
            :offset="[5, 5]">
            <ul>
              <li v-for="compound in wells[row][col]?.compounds" :key="compound.identifier">
                {{ compound.name }} ({{ compound.identifier }})
              </li>
            </ul>
          </q-tooltip>
        </td>
      </tr>
    </table>
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
</style>

<!-- :style="
            'background-color: ' + percentageToHsl(wells[row][col]?.measurements?[0]?.value || -1, 120, 0)
          " -->
