<script setup lang="ts">
import { computed, onMounted, ref, watch } from 'vue'
import ColorLegend from '~/components/plates/ColorLegend.vue'
import HeatMapSettings from '~/components/plates/HeatMapSettings.vue'
import PlateStats from '~/components/plates/PlateStats.vue'
import PlateTable from '~/components/plates/PlateTable.vue'
import { useExperimentStore } from '~/stores/experiments'
import { usePlateViewStore } from '~/stores/plateView'
import type { Plate } from '~/types/lab'
import { platePalettes } from '~/utils/plateHeatmap'
import { computeSSMD, computeZPrime } from '~/utils/plateStats'

type OverallStatsMap = {
  [key: string]: {
    min: number[]
    max: number[]
    mean: number[]
    median: number[]
    std: number[]
    mad: number[]
  }
}

const props = defineProps<{
  experimentId: number
  overallStats: OverallStatsMap
  availableMeasurementLabels: string[]
  timestamps: {
    [key: string]: string[]
  }
}>()

const emit = defineEmits<{
  (e: 'update-controls', payload: { pos: string | null; neg: string | null }): void
}>()

const { t } = useI18n()
const toast = useToast()
const experimentStore = useExperimentStore()
const plateViewStore = usePlateViewStore()

const experimentPlates = ref<Plate[]>([])
const isLoading = ref(false)

/**
 * Loads experiment plates:
 * - enables heatmap view,
 * - sorts plates by barcode,
 * - sets first measurement label by default.
 *
 * Loaded data example:
 * - `[{ id: 11, barcode: 'A001', ... }, { id: 12, barcode: 'A002', ... }]`
 */
const initialize = async () => {
  plateViewStore.showHeatmap = true
  plateViewStore.squareCompoundType = false
  plateViewStore.selectedMeasurement = props.availableMeasurementLabels[0] ?? null
  plateViewStore.selectedTimestampIdx = 0

  const savedControls = plateViewStore.getExperimentControls(props.experimentId)
  if (savedControls) {
    plateViewStore.selectedPosControl = savedControls.pos
    plateViewStore.selectedNegControl = savedControls.neg
  }

  isLoading.value = true
  try {
    const fetchedPlates = await experimentStore.fetchExperimentPlates(props.experimentId)
    const sortedPlates: Plate[] = []

    for (const plate of fetchedPlates) {
      sortedPlates.push(plate)
    }

    sortedPlates.sort((leftPlate, rightPlate) => leftPlate.barcode.localeCompare(rightPlate.barcode))
    experimentPlates.value = sortedPlates
  } finally {
    isLoading.value = false
  }
}

onMounted(async () => {
  await initialize()
})

watch(
  () => props.experimentId,
  async () => {
    await initialize()
  },
)

const timestampOptions = computed(() => {
  if (!plateViewStore.selectedMeasurement) {
    return []
  }

  const timestampValues = props.timestamps[plateViewStore.selectedMeasurement]
  if (!timestampValues || timestampValues.length === 0) {
    return []
  }

  const options: Array<{ label: string; value: number }> = []
  for (let index = 0; index < timestampValues.length; index += 1) {
    const timestamp = timestampValues[index]
    if (!timestamp) continue

    options.push({
      label: timestamp,
      value: index,
    })
  }

  return options
})

const max = computed(() => {
  if (!plateViewStore.selectedMeasurement) {
    return 0
  }

  const stats = props.overallStats[plateViewStore.selectedMeasurement]
  if (!stats || stats.max.length === 0) {
    return 0
  }

  return stats.max[0] ?? 0
})

const min = computed(() => {
  if (!plateViewStore.selectedMeasurement) {
    return 0
  }

  const stats = props.overallStats[plateViewStore.selectedMeasurement]
  if (!stats || stats.min.length === 0) {
    return 0
  }

  return stats.min[0] ?? 0
})

const getMinPerPlate = (plate: Plate) => {
  if (!plateViewStore.selectedMeasurement) {
    return 0
  }

  const measurementStats = plate.details.overall_stats[plateViewStore.selectedMeasurement]
  if (!measurementStats || measurementStats.min.length <= plateViewStore.selectedTimestampIdx) {
    return 0
  }

  return measurementStats.min[plateViewStore.selectedTimestampIdx] ?? 0
}

const getMaxPerPlate = (plate: Plate) => {
  if (!plateViewStore.selectedMeasurement) {
    return 0
  }

  const measurementStats = plate.details.overall_stats[plateViewStore.selectedMeasurement]
  if (!measurementStats || measurementStats.max.length <= plateViewStore.selectedTimestampIdx) {
    return 0
  }

  return measurementStats.max[plateViewStore.selectedTimestampIdx] ?? 0
}

const zPrimePerPlate = (plate: Plate) => {
  return computeZPrime(
    plate,
    plateViewStore.selectedMeasurement,
    plateViewStore.selectedPosControl,
    plateViewStore.selectedNegControl,
    plateViewStore.selectedTimestampIdx,
  )
}

const ssmdPerPlate = (plate: Plate) => {
  return computeSSMD(
    plate,
    plateViewStore.selectedMeasurement,
    plateViewStore.selectedPosControl,
    plateViewStore.selectedNegControl,
    plateViewStore.selectedTimestampIdx,
  )
}

const controlLabelOptions = computed(() => {
  const options: Array<{ label: string; value: string }> = []

  if (!plateViewStore.selectedMeasurement) {
    return options
  }

  if (experimentPlates.value.length === 0) {
    return options
  }

  const firstPlate = experimentPlates.value[0]
  if (!firstPlate) {
    return options
  }

  const statsForMeasurement = firstPlate.details.stats[plateViewStore.selectedMeasurement]
  if (!statsForMeasurement) {
    return options
  }

  const labelKeys = Object.keys(statsForMeasurement)
  for (const labelKey of labelKeys) {
    options.push({
      label: labelKey,
      value: labelKey,
    })
  }

  return options
})

const selectedPaletteLabel = computed({
  get: () => plateViewStore.heatmapPalette.label,
  set: (label: string) => {
    const palette = platePalettes.find((entry) => entry.label === label)
    if (!palette) return
    plateViewStore.heatmapPalette = palette
  },
})

const sendControlsToParent = () => {
  if (!plateViewStore.selectedPosControl || !plateViewStore.selectedNegControl) {
    toast.add({
      title: t('experiments.results.controls_missing'),
      color: 'error',
      duration: 3000,
    })
    return
  }

  plateViewStore.saveExperimentControls(
    props.experimentId,
    plateViewStore.selectedPosControl,
    plateViewStore.selectedNegControl,
  )

  emit('update-controls', {
    pos: plateViewStore.selectedPosControl,
    neg: plateViewStore.selectedNegControl,
  })

  toast.add({
    title: t('experiments.results.controls_saved'),
    color: 'success',
    duration: 2500,
  })
}

const openPlatePage = async (plate: Plate) => {
  await navigateTo(`/plates/${encodeURIComponent(plate.barcode)}`)
}

const canShowStats = (plate: Plate) => {
  const ssmdValue = ssmdPerPlate(plate)
  const zPrimeValue = zPrimePerPlate(plate)
  return ssmdValue !== null && zPrimeValue !== null
}
</script>

<template>
  <section>
    <HeatMapSettings :show-square-compound-type="false" :show-per-plate-view="true" />

    <p v-if="isLoading" class="text-sm text-slate-600">
      {{ t('experiments.results.loading') }}
    </p>
    <p v-else-if="experimentPlates.length === 0" class="text-sm text-slate-600">
      {{ t('experiments.results.no_plates') }}
    </p>

    <template v-else>
      <div v-if="props.availableMeasurementLabels.length > 0" class="mt-4 mb-6 flex flex-wrap items-end gap-3">
        <div class="w-[220px]">
          <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
            {{ t('plates.dynamic.controls.positive_control') }}
          </label>
          <select
            v-model="plateViewStore.selectedPosControl"
            class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-blue-300"
          >
            <option v-for="option in controlLabelOptions" :key="`pos-${option.value}`" :value="option.value">
              {{ option.label }}
            </option>
          </select>
        </div>

        <div class="w-[220px]">
          <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
            {{ t('plates.dynamic.controls.negative_control') }}
          </label>
          <select
            v-model="plateViewStore.selectedNegControl"
            class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-blue-300"
          >
            <option v-for="option in controlLabelOptions" :key="`neg-${option.value}`" :value="option.value">
              {{ option.label }}
            </option>
          </select>
        </div>

        <UButton
          color="secondary"
          variant="outline"
          :label="t('experiments.results.save_controls')"
          @click="sendControlsToParent"
        />
      </div>

      <div class="mt-2 mb-6 flex flex-wrap items-end gap-3">
        <div class="w-[240px]">
          <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
            {{ t('plates.dynamic.controls.select_label') }}
          </label>
          <select
            v-model="plateViewStore.selectedMeasurement"
            class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-blue-300"
          >
            <option v-for="label in props.availableMeasurementLabels" :key="`measurement-${label}`" :value="label">
              {{ label }}
            </option>
          </select>
        </div>

        <div v-if="timestampOptions.length > 0" class="w-[260px]">
          <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
            {{ t('plates.dynamic.controls.select_timestamp') }}
          </label>
          <select
            v-model.number="plateViewStore.selectedTimestampIdx"
            class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-blue-300"
          >
            <option v-for="option in timestampOptions" :key="`timestamp-${option.value}`" :value="option.value">
              {{ option.label }}
            </option>
          </select>
        </div>

        <div class="w-[220px]">
          <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
            {{ t('plates.dynamic.controls.select_color_palette') }}
          </label>
          <select
            v-model="selectedPaletteLabel"
            class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-blue-300"
          >
            <option v-for="palette in platePalettes" :key="`palette-${palette.label}`" :value="palette.label">
              {{ palette.label }}
            </option>
          </select>
        </div>
      </div>

      <div class="flex flex-wrap items-start justify-evenly gap-4">
        <article
          v-for="(plate, index) in experimentPlates"
          :key="`plate-${plate.id}`"
          class="cursor-pointer rounded-lg bg-white/70 p-3 shadow-sm"
          @click="openPlatePage(plate)"
        >
          <p class="mb-2 text-sm font-semibold text-blue-700">{{ plate.barcode }}</p>
          <div class="flex items-start gap-2">
            <PlateTable
              :plate-index="index"
              :plate="plate"
              :min="plateViewStore.perPlateView ? getMinPerPlate(plate) : min"
              :max="plateViewStore.perPlateView ? getMaxPerPlate(plate) : max"
            />
            <ColorLegend
              :min="plateViewStore.perPlateView ? getMinPerPlate(plate) : min"
              :max="plateViewStore.perPlateView ? getMaxPerPlate(plate) : max"
            />
          </div>

          <PlateStats v-if="canShowStats(plate)" :ssmd="ssmdPerPlate(plate)" :z-prime="zPrimePerPlate(plate)" />
        </article>
      </div>
    </template>
  </section>
</template>
