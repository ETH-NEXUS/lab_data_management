<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import BaseField from '~/components/common/BaseField.vue'
import WavesModalWrapper from '~/components/common/WavesModalWrapper.vue'
import { useExperimentStore } from '~/stores/experiments'
import type { BarcodeSide } from '~/types/experiments'
import { BARCODE_SIDE_OPTIONS } from '~/utils/barcode'
import { getErrorMessage } from '~/utils/errors'

const props = defineProps<{
  open: boolean
  experimentId: number
}>()

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
  (e: 'created'): void
}>()

const { t } = useI18n()
const toast = useToast()
const experimentStore = useExperimentStore()

const prefix = ref('')
const numberOfPlates = ref<number | null>(null)
const selectedSides = ref<BarcodeSide[]>([])

watch(
  () => props.open,
  (isOpen) => {
    if (!isOpen) return

    prefix.value = ''
    numberOfPlates.value = null
    selectedSides.value = []
  },
  { immediate: true },
)

/**
 * Maps numeric state to text field model and back.
 *
 * Data examples:
 * - state `{ numberOfPlates: 12 }` => field value `'12'`
 * - field value `''` => state `{ numberOfPlates: null }`
 */
const numberOfPlatesText = computed({
  get: (): string => {
    if (numberOfPlates.value === null) return ''
    return String(numberOfPlates.value)
  },
  set: (value: string) => {
    if (value.trim() === '') {
      numberOfPlates.value = null
      return
    }

    const parsedValue = Number(value)
    numberOfPlates.value = Number.isNaN(parsedValue) ? null : parsedValue
  },
})

/**
 * Validates form data before submit.
 *
 * Valid data example:
 * - `{ prefix: 'A001', numberOfPlates: 12, selectedSides: ['North', 'South'] }`
 */
const canSave = computed(() => {
  const hasPrefix = prefix.value.trim().length > 0
  const hasPlates = (numberOfPlates.value ?? 0) > 0
  const hasSides = selectedSides.value.length > 0
  return hasPrefix && hasPlates && hasSides && !experimentStore.isCreatingBarcodeSpecification
})

const close = () => emit('update:open', false)

/**
 * Custom `USelectMenu` UI to match the shared rounded field look used in modals.
 *
 * Input data example:
 * - `selectedSides = ['North', 'South']`
 *
 * Visual intent:
 * - rounded full trigger
 * - neutral border
 * - larger right chevron with extra spacing
 * - white dropdown surface
 */
const sidesSelectUi = {
  base: 'w-full rounded-full border border-slate-300 bg-white px-4 py-3 text-slate-600 outline-none ring-offset-0 transition focus-visible:ring-2 focus-visible:ring-blue-300',
  value: 'truncate text-slate-700',
  placeholder: 'truncate text-slate-400',
  trailing: 'pe-4',
  trailingIcon: 'size-6 text-slate-500',
  content: 'rounded-2xl bg-white shadow-md ring ring-slate-200',
  viewport: 'p-1',
  item: 'rounded-md',
}

const save = async () => {
  if (!canSave.value) return

  try {
    await experimentStore.createBarcodeSpecification({
      experiment_id: props.experimentId,
      prefix: prefix.value.trim(),
      number_of_plates: numberOfPlates.value ?? 0,
      sides: selectedSides.value,
    })

    emit('created')
    close()
    toast.add({
      title: t('experiments.barcodes.modal.created'),
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: t('experiments.barcodes.modal.create_failed'),
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  }
}
</script>

<template>
  <WavesModalWrapper
    :open="props.open"
    :title="t('experiments.barcodes.modal.title')"
    :description="t('experiments.barcodes.modal.description')"
    :dismissible="!experimentStore.isCreatingBarcodeSpecification"
    modal-class="w-full sm:max-w-3xl"
    body-container-class="w-full max-w-2xl px-8 pt-10"
    @update:open="emit('update:open', $event)"
  >
    <template #body>
      <div class="space-y-6">
        <BaseField
          v-model="prefix"
          :label="t('experiments.barcodes.fields.prefix')"
          :placeholder="t('experiments.barcodes.modal.prefix_placeholder')"
          :autofocus="true"
        />

        <BaseField
          v-model="numberOfPlatesText"
          :label="t('experiments.barcodes.fields.number_of_plates')"
          :placeholder="t('experiments.barcodes.modal.number_of_plates_placeholder')"
          type="number"
        />

        <div class="space-y-2">
          <label class="mb-1 block pl-4 text-sm font-medium">
            {{ t('experiments.barcodes.fields.sides') }}
          </label>
          <USelectMenu
            v-model="selectedSides"
            :items="BARCODE_SIDE_OPTIONS"
            value-key="value"
            label-key="label"
            multiple
            variant="none"
            class="w-full"
            :ui="sidesSelectUi"
            :placeholder="t('experiments.barcodes.modal.sides_placeholder')"
          />
        </div>
      </div>
    </template>

    <template #footer>
      <BaseButton
        :label="t('common.actions.cancel')"
        :on-click="close"
        variant="secondary"
        size="sm"
        width="auto"
        :disabled="experimentStore.isCreatingBarcodeSpecification"
      />
      <BaseButton
        :label="t('experiments.barcodes.modal.create_button')"
        :on-click="save"
        variant="primary"
        size="sm"
        width="auto"
        :loading="experimentStore.isCreatingBarcodeSpecification"
        :disabled="!canSave"
      />
    </template>
  </WavesModalWrapper>
</template>
