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

const canSubmit = computed(() => {
  if (Number.isNaN(parsedDmso.value)) return false
  if (Number.isNaN(parsedAmount.value)) return false
  if (parsedDmso.value < 0) return false
  if (parsedAmount.value < 0) return false
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
        <BaseField v-model="dmsoValue" :label="t('messages_page.modal.fields.dmso')" type="number" :autofocus="true" />

        <BaseField v-model="amountValue" :label="t('messages_page.modal.fields.amount')" type="number" />
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
