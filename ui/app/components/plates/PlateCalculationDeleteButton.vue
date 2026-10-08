<script setup lang="ts">
import { ref } from 'vue'
import { usePlateCalculations } from '~/composables/usePlateCalculations'
import type { Plate } from '~/types/lab'
import { getErrorMessage } from '~/utils/errors'

/**
 * Deletes a measurement calculated on the plate page, e.g. `Lum1_log10`, after
 * the user confirmed it. The server refuses an imported measurement.
 */
const props = defineProps<{
  plate: Plate
  label: string
}>()

const emit = defineEmits<{
  (e: 'deleted', label: string): void
}>()

const { t } = useI18n()
const toast = useToast()
const { isDeleting, deleteCalculation } = usePlateCalculations()

const isConfirmOpen = ref(false)

const openConfirmation = (): void => {
  isConfirmOpen.value = true
}

const closeConfirmation = (): void => {
  if (isDeleting.value) return
  isConfirmOpen.value = false
}

const onConfirmOpenChange = (isOpen: boolean) => {
  if (isDeleting.value) return
  isConfirmOpen.value = isOpen
}

const confirmDeletion = async () => {
  try {
    await deleteCalculation(props.plate.id, props.label)
    isConfirmOpen.value = false
    toast.add({ title: t('plates.calculations.delete.deleted', { label: props.label }), color: 'success' })
    emit('deleted', props.label)
  } catch (err: unknown) {
    toast.add({
      title: t('plates.calculations.delete.failed', { label: props.label }),
      description: getErrorMessage(err),
      color: 'error',
    })
  }
}
</script>

<template>
  <UButton
    color="error"
    variant="outline"
    size="sm"
    icon="i-heroicons-trash"
    :label="t('plates.calculations.delete.button')"
    @click="openConfirmation"
  />

  <!-- A plain confirmation card, as for archiving a plate -->
  <UModal
    :open="isConfirmOpen"
    :title="t('plates.calculations.delete.title', { label: props.label })"
    :description="t('plates.calculations.delete.description', { label: props.label, barcode: props.plate.barcode })"
    :dismissible="!isDeleting"
    class="w-full sm:max-w-lg"
    :ui="{ content: 'rounded-2xl bg-white shadow-md' }"
    @update:open="onConfirmOpenChange"
  >
    <template #footer>
      <div class="flex w-full justify-end gap-2">
        <UButton
          variant="ghost"
          color="neutral"
          :label="t('common.actions.cancel')"
          :disabled="isDeleting"
          @click="closeConfirmation"
        />
        <UButton
          color="error"
          :label="t('plates.calculations.delete.confirm')"
          :loading="isDeleting"
          :disabled="isDeleting"
          @click="confirmDeletion"
        />
      </div>
    </template>
  </UModal>
</template>
