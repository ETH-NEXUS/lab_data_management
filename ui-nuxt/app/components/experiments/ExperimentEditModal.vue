<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import { useExperimentStore } from '~/stores/experiments'
import type { ExperimentPayload } from '~/types/experiments'
import { getErrorMessage } from '~/utils/errors'

const props = defineProps<{
  open: boolean
  experimentId: number
  field: 'name' | 'description'
  initialValue: string
}>()

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
  (e: 'saved'): void
}>()

const toast = useToast()
const { t } = useI18n()
const experimentStore = useExperimentStore()
const localValue = ref('')

watch(
  () => [props.open, props.initialValue, props.field] as const,
  ([open]) => {
    if (open) {
      localValue.value = props.initialValue
    }
  },
  { immediate: true },
)

const title = computed(() => {
  return props.field === 'name' ? t('experiments.edit_modal.title_name') : t('experiments.edit_modal.title_description')
})

const description = computed(() => {
  return props.field === 'name'
    ? t('experiments.edit_modal.description_name')
    : t('experiments.edit_modal.description_description')
})

const canSave = computed(() => {
  if (props.field === 'name') {
    return localValue.value.trim().length > 0
  }

  return true
})

const close = () => emit('update:open', false)

/**
 * Sends a PATCH payload with one editable field.
 *
 * Payload examples:
 * - `{ name: 'Dose response v2' }`
 * - `{ description: 'Second run with adjusted controls' }`
 */
const save = async () => {
  if (!canSave.value) return

  const payload = { [props.field]: localValue.value } as ExperimentPayload

  try {
    await experimentStore.updateExperiment(props.experimentId, payload)
    emit('saved')
    close()
    toast.add({
      title: t('experiments.edit_modal.updated'),
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: t('experiments.edit_modal.update_failed'),
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  }
}
</script>

<template>
  <UModal
    :open="props.open"
    :title="title"
    :description="description"
    :dismissible="!experimentStore.isUpdatingExperiment"
    class="w-full sm:max-w-3xl"
    @update:open="emit('update:open', $event)"
  >
    <template #body>
      <div class="space-y-2">
        <p class="text-sm text-slate-600">
          {{
            props.field === 'name'
              ? t('experiments.edit_modal.field_name')
              : t('experiments.edit_modal.field_description')
          }}
        </p>

        <UInput
          v-if="props.field === 'name'"
          v-model="localValue"
          class="w-full"
          :placeholder="t('experiments.edit_modal.placeholder_name')"
          autofocus
        />
        <UTextarea
          v-else
          v-model="localValue"
          class="w-full"
          :placeholder="t('experiments.edit_modal.placeholder_description')"
          :rows="5"
          autofocus
        />
      </div>
    </template>

    <template #footer>
      <div class="flex w-full justify-end gap-2">
        <UButton color="neutral" variant="soft" :disabled="experimentStore.isUpdatingExperiment" @click="close">
          {{ t('common.actions.cancel') }}
        </UButton>
        <UButton
          color="primary"
          variant="solid"
          :loading="experimentStore.isUpdatingExperiment"
          :disabled="!canSave"
          @click="save"
        >
          {{ t('common.actions.save') }}
        </UButton>
      </div>
    </template>
  </UModal>
</template>
