<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import ColorLegend from '~/components/plates/ColorLegend.vue'
import PlateTable from '~/components/plates/PlateTable.vue'
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
const plateViewStore = usePlateViewStore()

// The corrected versions of the measurement shown in the main heatmap
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
</script>

<template>
  <section v-if="plateViewStore.showHeatmap && selectedCorrection" class="mt-8">
    <h3 class="mb-1 text-lg font-medium text-slate-800">
      {{ t('plates.background_correction.heatmap_title', { label: selectedCorrection.label }) }}
    </h3>
    <p class="mb-3 text-sm text-slate-600">
      {{
        t('plates.background_correction.formula', {
          method: t(`plates.background_correction.methods.${selectedCorrection.method}`),
          reference: selectedCorrection.referenceType,
        })
      }}
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

      <ColorLegend :min="minMax.min" :max="minMax.max" />
    </div>
  </section>
</template>
