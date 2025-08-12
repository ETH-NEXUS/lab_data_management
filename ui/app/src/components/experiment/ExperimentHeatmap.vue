<script setup lang="ts">
import {computed, onMounted, ref, PropType} from 'vue'
import {Plate} from '../models'
import {useProjectStore} from 'stores/project'
import PlateTable from 'components/plate/PlateTable.vue'
import HeatMapSettings from 'components/plate/HeatMapSettings.vue'
import ColorLegend from 'components/plate/ColorLegend.vue'
import {useI18n} from 'vue-i18n'
import {useRouter} from 'vue-router'
import {useQuasar} from 'quasar'
import {storeToRefs} from 'pinia'
import {useSettingsStore} from 'stores/settings'
import {palettes} from 'components/helpers'
import PlateStats from 'components/plate/PlateStats.vue'
import html2canvas from 'html2canvas'

const router = useRouter()
const $q = useQuasar()
const {platePage} = storeToRefs(useSettingsStore())

const selectedPosControl = ref<string | null>(null)
const selectedNegControl = ref<string | null>(null)

const emit = defineEmits(['updateControls'])

const sendControlsToParent = () => {
  if (!selectedPosControl.value || !selectedNegControl.value) {
    $q.notify({
      message: 'Please select both controls',
      color: 'negative',
      position: 'bottom',
    })
  }
  emit('updateControls', {
    pos: selectedPosControl.value,
    neg: selectedNegControl.value,
  })
}

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

const getMinPerPlate = (plate: Plate) => {
  if (selectedMeasurement.value) {
    return plate.details.overall_stats[selectedMeasurement.value].min[selectedTimestampIdx.value]
  }
  return 0
}

const getMaxPerPlate = (plate: Plate) => {
  if (selectedMeasurement.value) {
    return plate.details.overall_stats[selectedMeasurement.value].max[selectedTimestampIdx.value]
  }
  return 0
}

const z_primePerPlate = (plate: Plate) => {
  if (!selectedMeasurement.value) {
    return null
  }
  const measurement = selectedMeasurement.value
  const pos = selectedPosControl.value || 'P'
  const neg = selectedNegControl.value || 'N'

  if (!plate.details.stats) {
    return null
  }

  if (pos in plate.details.stats[measurement] && neg in plate.details.stats[measurement]) {
    const mad_pos = plate.details.stats[measurement][pos].mad[selectedTimestampIdx.value]
    const mad_neg = plate.details.stats[measurement][neg].mad[selectedTimestampIdx.value]
    const median_pos = plate.details.stats[measurement][pos].median[selectedTimestampIdx.value]
    const median_neg = plate.details.stats[measurement][neg].median[selectedTimestampIdx.value]

    return 1 - (3 * (mad_pos + mad_neg)) / Math.abs(median_pos - median_neg)
  }
  return null
}

const ssmdPerPlate = (plate: Plate) => {
  if (!selectedMeasurement.value) {
    return null
  }
  const measurement = selectedMeasurement.value
  const pos = selectedPosControl.value || 'P'
  const neg = selectedNegControl.value || 'N'

  if (!plate.details.stats) {
    return null
  }

  if (pos in plate.details.stats[measurement] && neg in plate.details.stats[measurement]) {
    const mad_pos = plate.details.stats[measurement][pos].mad[selectedTimestampIdx.value]
    const mad_neg = plate.details.stats[measurement][neg].mad[selectedTimestampIdx.value]
    const median_pos = plate.details.stats[measurement][pos].median[selectedTimestampIdx.value]
    const median_neg = plate.details.stats[measurement][neg].median[selectedTimestampIdx.value]

    return Math.abs(median_pos - median_neg) / (0.5 * (mad_pos + mad_neg))
  }
  return null
}

const controlLabelOptions = computed(() => {
  if (!selectedMeasurement.value) {
    return []
  }
  const statsForMeasurement = experimentPlates.value[0].details.stats[selectedMeasurement.value]
  return Object.keys(statsForMeasurement).map(labelKey => ({
    label: labelKey,
    value: labelKey,
  }))
})

const plateRefs = ref<HTMLElement[]>([])

const downloadPlateImage = async (index: number) => {
  const plateElement = plateRefs.value[index]
  if (!plateElement) return
  try {
    const canvas = await html2canvas(plateElement, {
      scale: 3, // makes it 3x higher resolution
    })
    const link = document.createElement('a')
    link.href = canvas.toDataURL('image/png')
    link.download = `plate-${index}.png`
    link.click()
  } catch (error) {
    console.error('Error generating image:', error)
  }
}
</script>

<template>
  <br />
  <HeatMapSettings :showSquareCompoundType="false" :show-per-plate-view="true" />
  <br />
  <div class="row q-mb-md q-my-xl" v-if="measurementOptions && measurementOptions.length > 0">
    <div class="col-3 q-mr-lg q-ml-lg">
      <q-select
        dense
        v-model="selectedPosControl"
        :options="controlLabelOptions"
        :label="t('label.positive_control')"
        emit-value
        map-options />
    </div>
    <div class="col-3">
      <q-select
        dense
        v-model="selectedNegControl"
        :options="controlLabelOptions"
        :label="t('label.negative_control')"
        emit-value
        map-options />
    </div>
    <q-btn class="q-ml-md" outline color="secondary" @click="sendControlsToParent">Save controls</q-btn>
  </div>

  <div class="col-4 q-mb-lg row">
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
      <div>
        <div ref="plateRefs" class="q-pa-md">
          <div class="q-mb-xs text-blue-8">{{ plate.barcode }}</div>
          <div class="fit row items-start content-start q-mr-lg">
            <PlateTable
              :plate-index="index"
              :plate="plate"
              :selectedMeasurement="selectedMeasurement"
              :selectedTimestampIdx="selectedTimestampIdx"
              :min="platePage.perPlateView ? getMinPerPlate(plate) : min"
              :max="platePage.perPlateView ? getMaxPerPlate(plate) : max" />
            <ColorLegend
              :max="platePage.perPlateView ? getMaxPerPlate(plate) : max"
              :min="platePage.perPlateView ? getMinPerPlate(plate) : min"
              :selectedMeasurement="selectedMeasurement" />
          </div>

          <PlateStats
            :ssmd="ssmdPerPlate(plate)"
            :z_prime="z_primePerPlate(plate)"
            v-if="ssmdPerPlate(plate) && z_primePerPlate(plate)" />
        </div>
        <q-btn color="primary" flat icon="download" @click.stop="downloadPlateImage(index)" class="q-ml-md" />
      </div>
    </div>
  </div>
</template>

<style scoped lang="sass"></style>
