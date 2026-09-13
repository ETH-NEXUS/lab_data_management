<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import BaseField from '~/components/common/BaseField.vue'
import WavesModalWrapper from '~/components/common/WavesModalWrapper.vue'
import type { Threshold, ThresholdUpdatePayload } from '~/types/messages'

const props = defineProps<{
  open: boolean
  threshold: Threshold
  isSaving?: boolean
  // Why the server refused the last save, shown below the fields.
  errorMessage?: string
}>()

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
  (e: 'submit', payload: ThresholdUpdatePayload): void
}>()

const { t } = useI18n()

const dmsoValue = ref('')
const amountValue = ref('')

/**
 * Syncs local form state when the modal is opened.
 *
 * Local data example:
 * - `{ dmsoValue: '80', amountValue: '2.5' }`
 */
watch(
  () => [props.open, props.threshold] as const,
  ([isOpen]) => {
    if (!isOpen) return

    dmsoValue.value = String(props.threshold.dmso)
    amountValue.value = String(props.threshold.amount)
  },
  { immediate: true },
)

const parsedDmso = computed(() => Number(dmsoValue.value))
const parsedAmount = computed(() => Number(amountValue.value))

/**
 * Explains what is wrong with the DMSO value, or returns null when it is fine.
 * An empty field counts as wrong: `Number('')` would otherwise save it as 0.
 *
 * Returned output examples:
 * - `'80'` -> `null`
 * - `'150'` -> `'DMSO must be a number between 0 and 100 %.'`
 */
const dmsoError = computed(() => {
  const isMissing = dmsoValue.value.trim() === '' || Number.isNaN(parsedDmso.value)
  if (isMissing || parsedDmso.value < 0 || parsedDmso.value > 100) {
    return t('messages_page.modal.errors.dmso_range')
  }
  return null
})

/**
 * Explains what is wrong with the volume value, or returns null when it is fine.
 *
 * Returned output examples:
 * - `'2.5'` -> `null`
 * - `'-1'` -> `'Amount must be a number of 0 µL or more.'`
 */
const amountError = computed(() => {
  const isMissing = amountValue.value.trim() === '' || Number.isNaN(parsedAmount.value)
  if (isMissing || parsedAmount.value < 0) {
    return t('messages_page.modal.errors.amount_min')
  }
  return null
})

const canSubmit = computed(() => {
  if (dmsoError.value || amountError.value) return false
  if (props.isSaving) return false
  return true
})

const close = () => emit('update:open', false)

/**
 * Emits validated threshold payload to the parent page.
 *
 * Payload example:
 * - `{ dmso: 75, amount: 2.0 }`
 */
const submit = () => {
  if (!canSubmit.value) return

  emit('submit', {
    dmso: parsedDmso.value,
    amount: parsedAmount.value,
  })
}
</script>

<template>
  <WavesModalWrapper
    :open="props.open"
    :title="t('messages_page.modal.title')"
    :description="t('messages_page.modal.description')"
    :dismissible="!props.isSaving"
    modal-class="w-full sm:max-w-3xl"
    body-container-class="w-full max-w-2xl px-8 pt-10"
    @update:open="emit('update:open', $event)"
  >
    <template #body>
      <div class="space-y-6">
        <div class="space-y-1">
          <BaseField
            v-model="dmsoValue"
            :label="t('messages_page.modal.fields.dmso')"
            type="number"
            :autofocus="true"
          />
          <p v-if="dmsoError" class="pl-1 text-xs text-red-600">{{ dmsoError }}</p>
        </div>

        <div class="space-y-1">
          <BaseField v-model="amountValue" :label="t('messages_page.modal.fields.amount')" type="number" />
          <p v-if="amountError" class="pl-1 text-xs text-red-600">{{ amountError }}</p>
        </div>

        <p v-if="props.errorMessage" class="text-sm text-red-600">{{ props.errorMessage }}</p>
      </div>
    </template>

    <template #footer>
      <BaseButton
        :label="t('common.actions.cancel')"
        :on-click="close"
        variant="secondary"
        size="sm"
        width="auto"
        :disabled="props.isSaving"
      />

      <BaseButton
        :label="t('common.actions.save')"
        :on-click="submit"
        variant="primary"
        size="sm"
        width="auto"
        :loading="props.isSaving"
        :disabled="!canSubmit"
      />
    </template>
  </WavesModalWrapper>
</template>
