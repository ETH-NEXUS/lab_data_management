<script setup lang="ts">
import {computed, onMounted, ref, PropType} from 'vue'
import {Plate} from '../models'
import {useProjectStore} from 'stores/project'
import PlateTable from 'components/plate/PlateTable.vue'
import HeatMapSettings from 'components/plate/HeatMapSettings.vue'
// import ColorLegend from 'components/ColorLegend.vue'
import {useI18n} from 'vue-i18n'
import {useRouter} from 'vue-router'
import {useQuasar} from 'quasar'
import {storeToRefs} from 'pinia'
import {useSettingsStore} from 'stores/settings'
import {palettes} from 'components/helpers'

const router = useRouter()
const $q = useQuasar()
const {platePage} = storeToRefs(useSettingsStore())

onMounted(async () => {
  if (platePage.value.showHeatmap == false) {
    platePage.value.showHeatmap = true
  }
  $q.loading.show({
    message: t('info.downloading_results'),
  })
  experimentPlates.value = await getExperimentPlates(props.experimentId)
  // sort by barcode
  experimentPlates.value.sort((a, b) => a.barcode.localeCompare(b.barcode))
  $q.loading.hide()
  measurementOptions.value = props.availableMeasurementLabels
  selectedMeasurement.value = measurementOptions.value[0]
})

const {getExperimentPlates} = useProjectStore()

const props = defineProps({
  experimentId: {
    type: Number,
    required: true,
  },
  overallStats: {
    type: Object as PropType<{
      [key: string]: {
        min: number[]
        max: number[]
        mean: number[]
        median: number[]
        std: number[]
        mad: number[]
      }
    }>,
    required: true,
  },
  availableMeasurementLabels: {
    type: Array as PropType<string[]>,
    required: true,
  },
  timestamps: {
    type: Object as PropType<{
      [key: string]: string[]
    }>,
    required: true,
  },
})
const {t} = useI18n()

const experimentPlates = ref<Plate[]>([])
const selectedMeasurement = ref<string | null>(null)
const measurementOptions = ref<string[]>([])
const selectedTimestampIdx = ref<number>(0)

const timestampOptions = computed(() => {
  if (selectedMeasurement.value) {
    const timestamps = props.timestamps[selectedMeasurement.value]
    if (timestamps && timestamps.length > 0) {
      return timestamps.map((ts, idx) => {
        return {label: ts, value: idx}
      })
    }
  }
  return null
})

const max = computed(() => {
  if (selectedMeasurement.value) {
    return props.overallStats[selectedMeasurement.value].max[0]
  }
  return 0
})

const min = computed(() => {
  if (selectedMeasurement.value) {
    return props.overallStats[selectedMeasurement.value].min[0]
  }
  return 0
})
</script>

<template>
  <br />
  <HeatMapSettings :showSquareCompoundType="false" />
  <br />

  <div class="col-4 q-mb-lg row">
    <!--    <div class="col-4">-->
    <!--      <ColorLegend :max="max" :min="min" :selectedMeasurement="selectedMeasurement" />-->
    <!--    </div>-->
    <div class="col-4 q-ml-lg">
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
      <q-select
        v-model="platePage.heatmapPalette"
        :options="palettes"
        :label="t('label.select_color_palette')"></q-select>
    </div>
  </div>
  <div class="fit row wrap justify-evenly items-start content-start">
    <div
      :style="{cursor: 'pointer'}"
      :key="`plate-${plate.id}`"
      v-for="(plate, index) in experimentPlates"
      class="q-mb-md q-ml-sm"
      @click="router.push(`/plate/${plate.barcode}`)">
      <div class="q-mb-xs text-blue-8">{{ plate.barcode }}</div>
      <PlateTable
        :plate-index="index"
        :plate="plate"
        :selectedMeasurement="selectedMeasurement"
        :selectedTimestampIdx="selectedTimestampIdx"
        :min="min"
        :max="max" />
    </div>
  </div>
</template>

<style scoped lang="sass"></style>
