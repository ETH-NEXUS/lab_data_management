<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import WavesModalWrapper from '~/components/common/WavesModalWrapper.vue'
import { usePlateBackgroundCorrection } from '~/composables/usePlateBackgroundCorrection'
import { usePlateLog10 } from '~/composables/usePlateLog10'
import { usePlateViewStore } from '~/stores/plateView'
import {
  BACKGROUND_CORRECTION_METHODS,
  PLATE_CALCULATIONS,
  type BackgroundCorrectionMethod,
  type PlateCalculation,
  type PlateCalculationResult,
} from '~/types/backgroundCorrection'
import type { Plate } from '~/types/lab'
import { getWellTypesOfMeasurement } from '~/utils/backgroundCorrection'
import { getErrorMessage } from '~/utils/errors'

const props = defineProps<{
  open: boolean
  plate: Plate
}>()

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
  (e: 'calculated', result: PlateCalculationResult): void
}>()

const { t } = useI18n()
const toast = useToast()
const plateViewStore = usePlateViewStore()
const { isCorrecting, correctPlateBackground } = usePlateBackgroundCorrection()
const { isCalculatingLog10, log10Measurement } = usePlateLog10()

const calculation = ref<PlateCalculation>('background_correction')
const label = ref<string | null>(null)
const referenceType = ref<string | null>(null)
const method = ref<BackgroundCorrectionMethod>('median')
const errorMessage = ref('')

const labels = computed(() => props.plate.details.measurement_labels ?? [])
const wellTypes = computed(() => getWellTypesOfMeasurement(props.plate, label.value))

// The reference wells the lab uses, in this order; "Nref" is imported as NREF
const USUAL_REFERENCE_TYPES = ['N1', 'NREF', 'R']

const defaultReferenceType = (): string | null => {
  const usual = USUAL_REFERENCE_TYPES.find((type) => wellTypes.value.includes(type))
  return usual ?? wellTypes.value[0] ?? null
}

watch(
  () => props.open,
  (isOpen) => {
    if (!isOpen) return

    calculation.value = 'background_correction'
    label.value = plateViewStore.selectedMeasurement ?? labels.value[0] ?? null
    referenceType.value = defaultReferenceType()
    method.value = 'median'
    errorMessage.value = ''
  },
  { immediate: true },
)

watch(label, () => {
  if (!referenceType.value || !wellTypes.value.includes(referenceType.value)) {
    referenceType.value = defaultReferenceType()
  }
})

const isRunning = computed(() => isCorrecting.value || isCalculatingLog10.value)

const canApply = computed(() => {
  if (isRunning.value || label.value === null) return false
  // Only the background correction needs reference wells
  return calculation.value === 'log10' || referenceType.value !== null
})

const close = () => emit('update:open', false)

const applyCorrection = async (sourceLabel: string, reference: string) => {
  const result = await correctPlateBackground(props.plate.id, {
    label: sourceLabel,
    reference_type: reference,
    method: method.value,
  })
  emit('calculated', { calculation: 'background_correction', label: sourceLabel, newLabel: result.label })
  toast.add({ title: t('plates.background_correction.success', { label: result.label }), color: 'success' })
}

const applyLog10 = async (sourceLabel: string) => {
  const result = await log10Measurement(props.plate.id, sourceLabel)
  emit('calculated', { calculation: 'log10', label: sourceLabel, newLabel: result.label })
  toast.add({
    title: t('plates.background_correction.success', { label: result.label }),
    description: result.skipped > 0 ? t('plates.calculations.log10_skipped', { count: result.skipped }) : undefined,
    color: 'success',
  })
}

const apply = async () => {
  if (!label.value) return

  errorMessage.value = ''
  try {
    if (calculation.value === 'log10') {
      await applyLog10(label.value)
    } else if (referenceType.value) {
      await applyCorrection(label.value, referenceType.value)
    }
    close()
  } catch (err: unknown) {
    errorMessage.value = getErrorMessage(err)
  }
}
</script>

<template>
  <WavesModalWrapper
    :open="props.open"
    :title="t('plates.calculations.title')"
    :description="t(`plates.calculations.descriptions.${calculation}`)"
    :dismissible="!isRunning"
    @update:open="emit('update:open', $event)"
  >
    <template #body>
      <div class="grid grid-cols-1 gap-4">
        <div>
          <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
            {{ t('plates.calculations.calculation') }}
          </label>
          <select
            v-model="calculation"
            class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-lime-500"
          >
            <option v-for="option in PLATE_CALCULATIONS" :key="`calculation-${option}`" :value="option">
              {{ t(`plates.calculations.names.${option}`) }}
            </option>
          </select>
        </div>

        <div>
          <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
            {{ t('plates.background_correction.measurement') }}
          </label>
          <select
            v-model="label"
            class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-lime-500"
          >
            <option v-for="option in labels" :key="`label-${option}`" :value="option">
              {{ option }}
            </option>
          </select>
        </div>

        <!-- log10 needs only the measurement -->
        <div v-if="calculation === 'background_correction'">
          <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
            {{ t('plates.background_correction.reference_type') }}
          </label>
          <select
            v-model="referenceType"
            class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-lime-500"
          >
            <option v-for="option in wellTypes" :key="`type-${option}`" :value="option">
              {{ option }}
            </option>
          </select>
        </div>

        <div v-if="calculation === 'background_correction'">
          <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
            {{ t('plates.background_correction.method') }}
          </label>
          <select
            v-model="method"
            class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-lime-500"
          >
            <option v-for="option in BACKGROUND_CORRECTION_METHODS" :key="`method-${option}`" :value="option">
              {{ t(`plates.background_correction.methods.${option}`) }}
            </option>
          </select>
        </div>

        <p v-if="errorMessage" class="rounded-lg border border-red-200 bg-red-50 p-3 text-sm text-red-700">
          {{ errorMessage }}
        </p>
      </div>
    </template>

    <template #footer>
      <BaseButton
        :label="t('common.actions.cancel')"
        :on-click="close"
        variant="secondary"
        size="sm"
        width="auto"
        :disabled="isRunning"
      />
      <BaseButton
        :label="t('plates.background_correction.apply_button')"
        :on-click="apply"
        variant="primary"
        size="sm"
        width="auto"
        :loading="isRunning"
        :disabled="!canApply"
      />
    </template>
  </WavesModalWrapper>
</template>
