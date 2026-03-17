<script setup lang="ts">
import { computed, watch } from 'vue'
import ColorLegend from '~/components/plates/ColorLegend.vue'
import HeatMapSettings from '~/components/plates/HeatMapSettings.vue'
import PlateStats from '~/components/plates/PlateStats.vue'
import PlateTable from '~/components/plates/PlateTable.vue'
import { usePlateViewStore } from '~/stores/plateView'
import type { Plate, WellInfo } from '~/types/lab'
import { platePalettes } from '~/utils/plateHeatmap'
import {
  computeSSMD,
  computeZPrime,
  getControlLabelOptions,
  getDefaultControlSelection,
  getOverallMinMaxForSelection,
  getTimestampOptions,
} from '~/utils/plateStats'

const props = defineProps<{
  plate: Plate
  minimalView?: boolean
}>()

const emit = defineEmits<{
  (e: 'well-selected', payload: WellInfo): void
  (e: 'refresh'): void
}>()

const { t } = useI18n()
const plateViewStore = usePlateViewStore()
const isMinimalView = computed(() => props.minimalView === true)

const measurementOptions = computed<string[]>(() => plateViewStore.measurementOptions ?? [])
const timestampOptions = computed(() => getTimestampOptions(props.plate, plateViewStore.selectedMeasurement))
const controlLabelOptions = computed(() => getControlLabelOptions(props.plate, plateViewStore.selectedMeasurement))

/**
 * Initializes measurement + control defaults when the opened plate changes.
 *
 * Behavior example:
 * - if labels are `['OD600', 'ATP']`, selects `OD600` initially.
 * - if stats contain `P1`/`N1`, picks them as control defaults.
 */
const initializeMeasurementState = (): void => {
  plateViewStore.measurementOptions = props.plate.details.measurement_labels ?? []
  plateViewStore.selectedMeasurement = measurementOptions.value[0] ?? null
  plateViewStore.selectedTimestampIdx = 0

  const defaults = getDefaultControlSelection(props.plate, plateViewStore.selectedMeasurement)
  plateViewStore.selectedPosControl = defaults.positive
  plateViewStore.selectedNegControl = defaults.negative
}

watch(
  () => props.plate.id,
  () => {
    initializeMeasurementState()
  },
  { immediate: true },
)

watch(
  () => plateViewStore.selectedMeasurement,
  () => {
    plateViewStore.selectedTimestampIdx = 0
    const defaults = getDefaultControlSelection(props.plate, plateViewStore.selectedMeasurement)
    plateViewStore.selectedPosControl = defaults.positive
    plateViewStore.selectedNegControl = defaults.negative
  },
)

const zPrime = computed(() =>
  computeZPrime(
    props.plate,
    plateViewStore.selectedMeasurement,
    plateViewStore.selectedPosControl,
    plateViewStore.selectedNegControl,
    plateViewStore.selectedTimestampIdx,
  ),
)

const ssmd = computed(() =>
  computeSSMD(
    props.plate,
    plateViewStore.selectedMeasurement,
    plateViewStore.selectedPosControl,
    plateViewStore.selectedNegControl,
    plateViewStore.selectedTimestampIdx,
  ),
)

const minMax = computed(() =>
  getOverallMinMaxForSelection(props.plate, plateViewStore.selectedMeasurement, plateViewStore.selectedTimestampIdx),
)

const wellContentOptions = computed(() => [
  { label: t('plates.dynamic.controls.hr_position'), value: 'hr_position' },
  { label: t('plates.dynamic.controls.well_type'), value: 'type' },
  { label: t('plates.dynamic.controls.index_position'), value: 'position' },
])

const selectedPaletteLabel = computed({
  get: () => plateViewStore.heatmapPalette.label,
  set: (label: string) => {
    const palette = platePalettes.find((entry) => entry.label === label)
    if (palette) {
      plateViewStore.heatmapPalette = palette
    }
  },
})

const onWellSelected = (wellInfo: WellInfo): void => {
  plateViewStore.selectedWellInfo = wellInfo
  emit('well-selected', wellInfo)
}
</script>

<template>
  <section>
    <HeatMapSettings v-if="!isMinimalView" />

    <div class="flex flex-nowrap gap-4">
      <div class="min-w-0 overflow-auto">
        <PlateTable :plate="props.plate" :min="minMax.min" :max="minMax.max" @well-selected="onWellSelected" />
      </div>

      <ColorLegend v-if="!isMinimalView" :min="minMax.min" :max="minMax.max" />
    </div>

    <PlateStats v-if="!isMinimalView && ssmd !== null && zPrime !== null" :ssmd="ssmd" :z-prime="zPrime" />

    <div v-if="!isMinimalView && measurementOptions.length > 0" class="mt-5 grid grid-cols-1 gap-4 lg:grid-cols-2">
      <div>
        <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
          {{ t('plates.dynamic.controls.positive_control') }}
        </label>
        <select
          v-model="plateViewStore.selectedPosControl"
          class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-lime-500"
        >
          <option v-for="option in controlLabelOptions" :key="`pos-${option.value}`" :value="option.value">
            {{ option.label }}
          </option>
        </select>
      </div>

      <div>
        <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
          {{ t('plates.dynamic.controls.negative_control') }}
        </label>
        <select
          v-model="plateViewStore.selectedNegControl"
          class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-lime-500"
        >
          <option v-for="option in controlLabelOptions" :key="`neg-${option.value}`" :value="option.value">
            {{ option.label }}
          </option>
        </select>
      </div>
    </div>

    <div v-if="isMinimalView || !props.plate.template" class="mt-6 grid grid-cols-1 gap-4 lg:grid-cols-2">
      <div>
        <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
          {{ t('plates.dynamic.controls.well_content') }}
        </label>
        <select
          v-model="plateViewStore.wellContent"
          class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-lime-500"
        >
          <option v-for="option in wellContentOptions" :key="`well-content-${option.value}`" :value="option.value">
            {{ option.label }}
          </option>
        </select>
      </div>

      <div v-if="!isMinimalView && measurementOptions.length > 0">
        <label class="inline-flex cursor-pointer items-center gap-2 text-sm text-slate-700">
          <input v-model="plateViewStore.showHeatmap" type="checkbox" class="cursor-pointer" />
          <span>{{ t('plates.dynamic.controls.show_heatmap') }}</span>
        </label>

        <div v-if="plateViewStore.showHeatmap" class="mt-3 grid grid-cols-1 gap-3">
          <div>
            <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
              {{ t('plates.dynamic.controls.select_label') }}
            </label>
            <select
              v-model="plateViewStore.selectedMeasurement"
              class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-lime-500"
            >
              <option v-for="label in measurementOptions" :key="`measurement-${label}`" :value="label">
                {{ label }}
              </option>
            </select>
          </div>

          <div v-if="timestampOptions.length > 0">
            <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
              {{ t('plates.dynamic.controls.select_timestamp') }}
            </label>
            <select
              v-model.number="plateViewStore.selectedTimestampIdx"
              class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-lime-500"
            >
              <option v-for="option in timestampOptions" :key="`timestamp-${option.value}`" :value="option.value">
                {{ option.label }}
              </option>
            </select>
          </div>

          <div>
            <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
              {{ t('plates.dynamic.controls.select_color_palette') }}
            </label>
            <select
              v-model="selectedPaletteLabel"
              class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-lime-500"
            >
              <option v-for="palette in platePalettes" :key="`palette-${palette.label}`" :value="palette.label">
                {{ palette.label }}
              </option>
            </select>
          </div>
        </div>
      </div>
    </div>
  </section>
</template>
