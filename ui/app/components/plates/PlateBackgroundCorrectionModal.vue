<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import WavesModalWrapper from '~/components/common/WavesModalWrapper.vue'
import { usePlateBackgroundCorrection } from '~/composables/usePlateBackgroundCorrection'
import { usePlateViewStore } from '~/stores/plateView'
import { BACKGROUND_CORRECTION_METHODS, type BackgroundCorrectionMethod } from '~/types/backgroundCorrection'
import type { Plate } from '~/types/lab'
import { getWellTypesOfMeasurement } from '~/utils/backgroundCorrection'
import { getErrorMessage } from '~/utils/errors'

const props = defineProps<{
  open: boolean
  plate: Plate
}>()

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
  // The measurement that was corrected and the new one, e.g. 'Lum_CTG', 'Lum_CTG_bc_N1_median'
  (e: 'corrected', label: string, correctedLabel: string): void
}>()

const { t } = useI18n()
const toast = useToast()
const plateViewStore = usePlateViewStore()
const { isCorrecting, correctPlateBackground } = usePlateBackgroundCorrection()

const label = ref<string | null>(null)
const referenceType = ref<string | null>(null)
const method = ref<BackgroundCorrectionMethod>('median')
const errorMessage = ref('')

const labels = computed(() => props.plate.details.measurement_labels ?? [])
const wellTypes = computed(() => getWellTypesOfMeasurement(props.plate, label.value))

// The reference wells the lab uses, in this order; "Nref" is imported as NREF
const USUAL_REFERENCE_TYPES = ['N1', 'NREF']

const defaultReferenceType = (): string | null => {
  const usual = USUAL_REFERENCE_TYPES.find((type) => wellTypes.value.includes(type))
  return usual ?? wellTypes.value[0] ?? null
}

watch(
  () => props.open,
  (isOpen) => {
    if (!isOpen) return

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

const canApply = computed(() => {
  return !isCorrecting.value && label.value !== null && referenceType.value !== null
})

const close = () => emit('update:open', false)

const applyCorrection = async () => {
  if (!label.value || !referenceType.value) return

  errorMessage.value = ''
  try {
    const result = await correctPlateBackground(props.plate.id, {
      label: label.value,
      reference_type: referenceType.value,
      method: method.value,
    })
    emit('corrected', label.value, result.label)
    close()
    toast.add({
      title: t('plates.background_correction.success', { label: result.label }),
      color: 'success',
      duration: 3000,
    })
  } catch (err: unknown) {
    errorMessage.value = getErrorMessage(err)
  }
}
</script>

<template>
  <WavesModalWrapper
    :open="props.open"
    :title="t('plates.background_correction.title')"
    :description="t('plates.background_correction.description')"
    :dismissible="!isCorrecting"
    @update:open="emit('update:open', $event)"
  >
    <template #body>
      <div class="grid grid-cols-1 gap-4">
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

        <div>
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

        <div>
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
        :disabled="isCorrecting"
      />
      <BaseButton
        :label="t('plates.background_correction.apply_button')"
        :on-click="applyCorrection"
        variant="primary"
        size="sm"
        width="auto"
        :loading="isCorrecting"
        :disabled="!canApply"
      />
    </template>
  </WavesModalWrapper>
</template>
