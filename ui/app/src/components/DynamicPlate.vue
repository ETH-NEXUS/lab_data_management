<script setup lang="ts">
import {computed, defineProps, defineEmits, PropType, ref, onMounted} from 'vue'
import {Plate, Well} from './models'
import {positionFromRowCol} from '../helpers/plate'
import {palettes} from 'components/data'
import {percentageToHsl} from 'components/helpers'

const props = defineProps({
  plate: {
    type: Object as PropType<Plate>,
    required: true,
  },
})

onMounted(() => {
  measurementsValuesIndices.value = findNumberOfMeasurements()
  measurementsMetadata.value = findMeasurementMetadata(currentMeasurementValueIndex.value)
})

const currentMeasurementValueIndex = ref<number>(0)
const showHeatmap = ref<boolean>(false)
const measurementsValuesIndices = ref<number[]>([])
const measurementsMetadata = ref<{feature: string; abbreviation: string; unit: string} | null>(null)
const palette = ref<string>('green_red')

type PaletteToFunction = {
  [key in keyof typeof palettes]: (percentage: number) => string
}

const paletteToFunction: PaletteToFunction = {} as PaletteToFunction

for (const key in palettes) {
  if (palettes.hasOwnProperty(key)) {
    const {from, to} = palettes[key]
    paletteToFunction[key] = (percentage: number) => percentageToHsl(percentage, from, to)
  }
}

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

const maxMeasurement = computed(() => {
  return findMaxMeasurement()
})

const minMeasurement = computed(() => {
  return findMinMeasurement()
})

const findNumberOfMeasurements = () => {
  const {wells = []} = props.plate
  const maxMeasurements = Math.max(...wells.map(({measurements = []}) => measurements.length), 0)
  return Array.from({length: maxMeasurements}, (_, index) => index)
}

const findMaxMeasurement = () => {
  return props.plate.wells?.reduce((max, {measurements = []}) => {
    return Math.max(max, ...measurements.map(m => m.value))
  }, 0)
}
const findMinMeasurement = () => {
  return props.plate.wells?.reduce((min, {measurements = []}) => {
    return Math.min(min, ...measurements.map(m => m.value))
  }, 0)
}

const findMeasurementPercentage = (value: number) => {
  if (minMeasurement.value !== undefined && maxMeasurement.value !== undefined && value) {
    return (value - minMeasurement.value) / (maxMeasurement.value - minMeasurement.value)
  }
  return 0
}

const findMeasurementMetadata = (n: number) => {
  const well = props.plate?.wells?.[0]
  const measurement = well?.measurements?.[n]
  if (measurement) {
    const {name, abbrev, unit} = measurement.feature
    return {
      feature: name === null ? 'not specified' : name,
      abbreviation: abbrev === null ? 'not specified' : abbrev,
      unit: unit === null ? 'not specified' : unit,
    }
  }
  return null
}

const toggleHeatmap = () => {
  showHeatmap.value = !showHeatmap.value
}

const changeCurrentValueIndex = (n: number) => {
  currentMeasurementValueIndex.value = n
}
</script>

<template>
  <div v-if="measurementsValuesIndices.length > 0" class="q-pa-sm">
    <q-btn flat color="primary" @click="toggleHeatmap">Show Measurement Heatmap</q-btn>
    <div class="q-pa-sm" v-if="measurementsValuesIndices.length > 1">
      <p
        @click="changeCurrentValueIndex(n)"
        class="text-blue-5 cursor-pointer"
        v-for="n of measurementsValuesIndices"
        :key="`index_${n}`">
        >> {{ `Value ${n}` }}
      </p>
    </div>
    <div v-if="measurementsMetadata && showHeatmap" class="q-ml-md q-mt-sm">
      <p>
        <b>Feature</b>
        : {{ measurementsMetadata.feature }}
      </p>
      <p>
        <b>Abbreviation</b>
        : {{ measurementsMetadata.abbreviation }}
      </p>
      <p>
        <b>Unit</b>
        : {{ measurementsMetadata.unit }}
      </p>
    </div>

    <div class="q-my-md" v-if="showHeatmap">
      <q-radio
        :key="p.label"
        v-for="p of Object.values(palettes)"
        v-model="palette"
        :val="p.val"
        :label="p.label"></q-radio>
    </div>
  </div>

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
          :style="{
            backgroundColor: showHeatmap
              ? paletteToFunction[palette](
                  findMeasurementPercentage(
                    wells[row][col]?.measurements[currentMeasurementValueIndex]?.value
                      ? wells[row][col]?.measurements[currentMeasurementValueIndex].value
                      : 0
                  ),
                  120,
                  0
                )
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
