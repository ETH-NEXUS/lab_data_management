<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import ColorLegend from '~/components/plates/ColorLegend.vue'
import PlateBackgroundCorrectionModal from '~/components/plates/PlateBackgroundCorrectionModal.vue'
import PlateTable from '~/components/plates/PlateTable.vue'
import { usePlateStore } from '~/stores/plates'
import { usePlateViewStore } from '~/stores/plateView'
import type { Plate, WellInfo } from '~/types/lab'
import { findBackgroundCorrections } from '~/utils/backgroundCorrection'
import { getOverallMinMaxForSelection } from '~/utils/plateStats'

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

// The corrected versions of the measurement chosen for the main heatmap
const corrections = computed(() => findBackgroundCorrections(props.plate, plateViewStore.selectedMeasurement))

const selectedLabel = ref<string | null>(null)

watch(
  corrections,
  () => {
    const labels = corrections.value.map((correction) => correction.label)
    if (!selectedLabel.value || !labels.includes(selectedLabel.value)) {
      selectedLabel.value = labels[0] ?? null
    }
  },
  { immediate: true },
)

const selectedCorrection = computed(() => {
  return corrections.value.find((correction) => correction.label === selectedLabel.value) ?? null
})

// The corrected measurement has the same time points as the original one
const minMax = computed(() =>
  getOverallMinMaxForSelection(props.plate, selectedLabel.value, plateViewStore.selectedTimestampIdx),
)

/**
 * Reloads only the plate (not the whole page), so the heatmap settings stay,
 * and shows the corrected measurement: the main heatmap gets its original one.
 *
 * Accepted input example: `'Lum_CTG'`
 */
const onCorrected = async (label: string): Promise<void> => {
  try {
    const plate = await plateStore.fetchPlateByBarcode(props.plate.barcode)
    if (!plate) return

    plateViewStore.measurementOptions = plate.details.measurement_labels ?? []
    plateViewStore.selectedMeasurement = label
  } catch (err) {
    // The plate store keeps the error and the page shows it
    console.error(err)
  }
}
</script>

<template>
  <section class="mt-8 rounded-xl border border-black/10 bg-white/50 p-4">
    <div class="flex flex-wrap items-center justify-between gap-3">
      <h3 class="text-lg font-medium text-slate-800">{{ t('plates.background_correction.section_title') }}</h3>
      <UButton
        color="secondary"
        variant="outline"
        icon="i-heroicons-adjustments-horizontal"
        :label="t('plates.background_correction.open_button')"
        @click="isModalOpen = true"
      />
    </div>

    <p v-if="!selectedCorrection" class="mt-2 text-sm text-slate-600">
      {{ t('plates.background_correction.none_yet', { label: plateViewStore.selectedMeasurement ?? '' }) }}
    </p>

    <template v-else>
      <h4 class="mt-3 mb-1 font-medium text-slate-800">
        {{ t('plates.background_correction.heatmap_title', { label: selectedCorrection.label }) }}
      </h4>
      <p class="text-sm text-slate-600">{{ t('plates.background_correction.formula_caption') }}</p>
      <!-- Monospace with spaces, so "-" reads as a minus and not as a dash -->
      <code class="my-1 block w-fit rounded-md bg-slate-100 px-3 py-2 font-mono text-sm text-slate-800">
        {{
          t('plates.background_correction.formula', {
            method: t(`plates.background_correction.methods.${selectedCorrection.method}`),
            reference: selectedCorrection.referenceType,
          })
        }}
      </code>
      <p class="mb-3 text-xs text-slate-500">
        {{ t('plates.background_correction.formula_note', { reference: selectedCorrection.referenceType }) }}
      </p>

      <div v-if="corrections.length > 1" class="mb-3 max-w-sm">
        <select
          v-model="selectedLabel"
          class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-lime-500"
        >
          <option v-for="correction in corrections" :key="`correction-${correction.label}`" :value="correction.label">
            {{ correction.label }}
          </option>
        </select>
      </div>

      <!-- Always a heatmap, also without "Show heatmap" for the main one -->
      <div class="flex flex-nowrap gap-4">
        <div class="min-w-0 overflow-auto">
          <PlateTable
            :plate="props.plate"
            :min="minMax.min"
            :max="minMax.max"
            :measurement-label="selectedCorrection.label"
            @well-selected="emit('well-selected', $event)"
          />
        </div>

        <ColorLegend :min="minMax.min" :max="minMax.max" always-shown />
      </div>
    </template>

    <PlateBackgroundCorrectionModal v-model:open="isModalOpen" :plate="props.plate" @corrected="onCorrected" />
  </section>
</template>
