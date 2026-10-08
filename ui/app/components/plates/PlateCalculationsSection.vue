<script setup lang="ts">
import { computed, ref } from 'vue'
import PlateCalculationCard from '~/components/plates/PlateCalculationCard.vue'
import PlateCalculationModal from '~/components/plates/PlateCalculationModal.vue'
import PlateRawDataExplanation from '~/components/plates/PlateRawDataExplanation.vue'
import { usePlateStore } from '~/stores/plates'
import { usePlateViewStore } from '~/stores/plateView'
import type { Plate, WellInfo } from '~/types/lab'
import { getCalculatedLabels } from '~/utils/plateDatasets'

/**
 * The "Calculations" of the plate page: the button to start one, what the main
 * heatmap above shows, and every measurement calculated on this page (log10,
 * background corrections, %Inhibition, %Activity) with its own heatmap, its
 * explanation and a button to delete it.
 */
const props = defineProps<{
  plate: Plate
}>()

const emit = defineEmits<{
  (e: 'well-selected', payload: WellInfo): void
}>()

const { t } = useI18n()
const plateStore = usePlateStore()
const plateViewStore = usePlateViewStore()

const isModalOpen = ref(false)

const openModal = (): void => {
  isModalOpen.value = true
}

// e.g. ['Lum1_log10', 'Lum1_bc_R_median']
const calculatedLabels = computed(() => getCalculatedLabels(props.plate))

// The main heatmap shows a calculated measurement: it is explained with its own heatmap below
const isMainHeatmapCalculated = computed(() => {
  const label = plateViewStore.selectedMeasurement
  return label !== null && calculatedLabels.value.includes(label)
})

/**
 * Reloads only the plate (not the whole page) after a calculation or a deletion,
 * so the heatmap settings stay. If the main heatmap showed a deleted measurement,
 * it shows the first one of the plate instead.
 */
const reloadPlate = async (): Promise<void> => {
  try {
    const plate = await plateStore.fetchPlateByBarcode(props.plate.barcode)
    if (!plate) return

    const labels = plate.details.measurement_labels ?? []
    plateViewStore.measurementOptions = labels
    if (plateViewStore.selectedMeasurement && !labels.includes(plateViewStore.selectedMeasurement)) {
      plateViewStore.selectedMeasurement = labels[0] ?? null
    }
    // The well details on the right show the changed measurements too
    const selectedPosition = plateViewStore.selectedWellInfo?.position
    if (selectedPosition !== undefined) {
      const selectedWell = plate.wells.find((well) => well.position === selectedPosition)
      plateViewStore.selectedWellInfo = { well: selectedWell, position: selectedPosition }
    }
  } catch (err) {
    // The plate store keeps the error and the page shows it
    console.error(err)
  }
}
</script>

<template>
  <section class="mt-8 rounded-xl border border-black/10 bg-white/50 p-4">
    <div class="flex flex-wrap items-center justify-between gap-3">
      <h3 class="text-lg font-medium text-slate-800">{{ t('plates.calculations.section_title') }}</h3>
      <UButton
        color="secondary"
        variant="outline"
        icon="i-heroicons-adjustments-horizontal"
        :label="t('plates.calculations.open_button')"
        @click="openModal"
      />
    </div>

    <p v-if="isMainHeatmapCalculated" class="mt-2 text-sm text-slate-700">
      {{ t('plates.calculations.main_is_calculated', { label: plateViewStore.selectedMeasurement }) }}
    </p>
    <PlateRawDataExplanation
      v-else-if="plateViewStore.selectedMeasurement"
      :plate="props.plate"
      :label="plateViewStore.selectedMeasurement"
      :timestamp-index="plateViewStore.selectedTimestampIdx"
    />

    <p v-if="calculatedLabels.length === 0" class="mt-4 text-sm text-slate-600">
      {{ t('plates.calculations.none_calculated') }}
    </p>
    <PlateCalculationCard
      v-for="label in calculatedLabels"
      :key="`calculation-${label}`"
      :plate="props.plate"
      :label="label"
      @well-selected="emit('well-selected', $event)"
      @deleted="reloadPlate"
    />

    <PlateCalculationModal v-model:open="isModalOpen" :plate="props.plate" @calculated="reloadPlate" />
  </section>
</template>
