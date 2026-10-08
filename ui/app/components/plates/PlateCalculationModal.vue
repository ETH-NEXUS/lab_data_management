<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import WavesModalWrapper from '~/components/common/WavesModalWrapper.vue'
import { usePlateCalculations } from '~/composables/usePlateCalculations'
import { usePlateViewStore } from '~/stores/plateView'
import type { Plate } from '~/types/lab'
import {
  BACKGROUND_CORRECTION_METHODS,
  PLATE_CALCULATIONS,
  type BackgroundCorrectionMethod,
  type PlateCalculation,
  type PlateCalculationResult,
  type PlateCalculationSettings,
} from '~/types/plateCalculations'
import { getErrorMessage } from '~/utils/errors'
import { getWellTypesOfMeasurement } from '~/utils/plateDatasets'

/**
 * The window to start a calculation with one measurement of the plate:
 * background correction, log10 or the normalization (%Inhibition, %Activity).
 */
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
const { isCalculating, calculate } = usePlateCalculations()

const SELECT_CLASS =
  'w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-lime-500'

// The well types the lab uses for each role, in this order; "Nref" is imported as NREF
const USUAL_REFERENCE_TYPES = ['N1', 'NREF', 'R']
const USUAL_NEGATIVE_TYPES = ['N', 'N1']
const USUAL_POSITIVE_TYPES = ['P', 'P1']

const calculation = ref<PlateCalculation>('background_correction')
const label = ref<string | null>(null)
const referenceType = ref<string | null>(null)
const method = ref<BackgroundCorrectionMethod>('median')
const negativeType = ref<string | null>(null)
const positiveType = ref<string | null>(null)
const errorMessage = ref('')

const labels = computed(() => props.plate.details.measurement_labels ?? [])
const wellTypes = computed(() => getWellTypesOfMeasurement(props.plate, label.value))

const usualWellType = (usualTypes: string[]): string | null => {
  return usualTypes.find((type) => wellTypes.value.includes(type)) ?? null
}

// A chosen well type that the measurement does not have is chosen again
const keepOrChooseWellType = (chosen: string | null, usualTypes: string[], fallback: string | null) => {
  if (chosen && wellTypes.value.includes(chosen)) return chosen
  return usualWellType(usualTypes) ?? fallback
}

const chooseWellTypes = (): void => {
  referenceType.value = keepOrChooseWellType(referenceType.value, USUAL_REFERENCE_TYPES, wellTypes.value[0] ?? null)
  negativeType.value = keepOrChooseWellType(negativeType.value, USUAL_NEGATIVE_TYPES, null)
  positiveType.value = keepOrChooseWellType(positiveType.value, USUAL_POSITIVE_TYPES, null)
}

watch(
  () => props.open,
  (isOpen) => {
    if (!isOpen) return

    calculation.value = 'background_correction'
    label.value = plateViewStore.selectedMeasurement ?? labels.value[0] ?? null
    referenceType.value = null
    negativeType.value = null
    positiveType.value = null
    chooseWellTypes()
    method.value = 'median'
    errorMessage.value = ''
  },
  { immediate: true },
)

watch(label, chooseWellTypes)

/**
 * What is sent to the server, or null while something is missing.
 *
 * Returned data example (log10):
 * - `{ label: 'Lum1' }`
 */
const settings = computed((): PlateCalculationSettings | null => {
  if (!label.value) return null
  if (calculation.value === 'log10') return { label: label.value }
  if (calculation.value === 'normalization') {
    if (!negativeType.value || !positiveType.value || negativeType.value === positiveType.value) return null
    return { label: label.value, negative_type: negativeType.value, positive_type: positiveType.value }
  }
  if (!referenceType.value) return null
  return { label: label.value, reference_type: referenceType.value, method: method.value }
})

const canApply = computed(() => !isCalculating.value && settings.value !== null)

const close = () => emit('update:open', false)

const apply = async () => {
  if (!settings.value || !label.value) return

  errorMessage.value = ''
  try {
    const result = await calculate(props.plate.id, calculation.value, settings.value)
    emit('calculated', { calculation: calculation.value, label: label.value, newLabel: result.label })
    close()
    // The normalization saves two measurements: the %Inhibition and the %Activity
    const title = result.activity_label
      ? t('plates.calculations.success_two', { label: result.label, activityLabel: result.activity_label })
      : t('plates.calculations.success', { label: result.label })
    toast.add({
      title,
      // Only log10 and the normalization leave wells out
      description: result.skipped ? t('plates.calculations.log10.skipped', { count: result.skipped }) : undefined,
      color: 'success',
    })
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
    :dismissible="!isCalculating"
    @update:open="emit('update:open', $event)"
  >
    <template #body>
      <div class="grid grid-cols-1 gap-4">
        <div>
          <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
            {{ t('plates.calculations.calculation') }}
          </label>
          <select v-model="calculation" :class="SELECT_CLASS">
            <option v-for="option in PLATE_CALCULATIONS" :key="`calculation-${option}`" :value="option">
              {{ t(`plates.calculations.names.${option}`) }}
            </option>
          </select>
        </div>

        <div>
          <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
            {{ t('plates.calculations.measurement') }}
          </label>
          <select v-model="label" :class="SELECT_CLASS">
            <option v-for="option in labels" :key="`label-${option}`" :value="option">
              {{ option }}
            </option>
          </select>
        </div>

        <template v-if="calculation === 'background_correction'">
          <div>
            <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
              {{ t('plates.background_correction.reference_type') }}
            </label>
            <select v-model="referenceType" :class="SELECT_CLASS">
              <option v-for="option in wellTypes" :key="`reference-${option}`" :value="option">
                {{ option }}
              </option>
            </select>
          </div>

          <div>
            <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
              {{ t('plates.background_correction.method') }}
            </label>
            <select v-model="method" :class="SELECT_CLASS">
              <option v-for="option in BACKGROUND_CORRECTION_METHODS" :key="`method-${option}`" :value="option">
                {{ t(`plates.background_correction.methods.${option}`) }}
              </option>
            </select>
          </div>
        </template>

        <template v-if="calculation === 'normalization'">
          <div>
            <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
              {{ t('plates.calculations.negative_control') }}
            </label>
            <select v-model="negativeType" :class="SELECT_CLASS">
              <option v-for="option in wellTypes" :key="`negative-${option}`" :value="option">
                {{ option }}
              </option>
            </select>
          </div>

          <div>
            <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
              {{ t('plates.calculations.positive_control') }}
            </label>
            <select v-model="positiveType" :class="SELECT_CLASS">
              <option v-for="option in wellTypes" :key="`positive-${option}`" :value="option">
                {{ option }}
              </option>
            </select>
          </div>
        </template>

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
        :disabled="isCalculating"
      />
      <BaseButton
        :label="t('plates.calculations.apply_button')"
        :on-click="apply"
        variant="primary"
        size="sm"
        width="auto"
        :loading="isCalculating"
        :disabled="!canApply"
      />
    </template>
  </WavesModalWrapper>
</template>
