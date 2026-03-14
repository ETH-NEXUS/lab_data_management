<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import { useExperimentStore } from '~/stores/experiments'
import type { Experiment } from '~/types/experiments'
import { getErrorMessage } from '~/utils/errors'

const props = defineProps<{
  open: boolean
  projectId: number | null
  projectName?: string | null
}>()

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
  (e: 'created', experiment: Experiment): void
}>()

const toast = useToast()
const { t } = useI18n()
const experimentStore = useExperimentStore()
const experimentName = ref('')

watch(
  () => props.open,
  (isOpen) => {
    if (isOpen) {
      experimentName.value = ''
    }
  },
  { immediate: true },
)

const modalTitle = computed(() => {
  if (props.projectName) {
    return t('experiments.create_modal.title_with_project', { projectName: props.projectName })
  }

  return t('experiments.create_modal.title')
})

const canSave = computed(() => {
  return experimentName.value.trim().length > 0 && (props.projectId ?? 0) > 0 && !experimentStore.isCreatingExperiment
})

const close = () => emit('update:open', false)

/**
 * Creates a new experiment for the selected project.
 *
 * Payload example sent by store method:
 * - `{ name: 'Dose response', project: 5 }`
 */
const save = async () => {
  if (!canSave.value || !props.projectId) return

  try {
    const experiment = await experimentStore.createExperiment({
      name: experimentName.value.trim(),
      project: props.projectId,
    })

    emit('created', experiment)
    close()
    toast.add({
      title: t('experiments.create_modal.created'),
      description: experiment.name,
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: t('experiments.create_modal.create_failed'),
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
    :title="modalTitle"
    :description="t('experiments.create_modal.description')"
    :dismissible="!experimentStore.isCreatingExperiment"
    class="w-full sm:max-w-2xl"
    @update:open="emit('update:open', $event)"
  >
    <template #body>
      <div class="space-y-2">
        <p class="text-sm text-slate-600">{{ t('experiments.create_modal.field_name') }}</p>
        <UInput
          v-model="experimentName"
          class="w-full"
          :placeholder="t('experiments.create_modal.placeholder_name')"
          autofocus
        />
      </div>
    </template>

    <template #footer>
      <div class="flex w-full justify-end gap-2">
        <UButton color="neutral" variant="soft" :disabled="experimentStore.isCreatingExperiment" @click="close">
          {{ t('common.actions.cancel') }}
        </UButton>
        <UButton
          color="primary"
          variant="solid"
          :loading="experimentStore.isCreatingExperiment"
          :disabled="!canSave"
          @click="save"
        >
          {{ t('experiments.create_modal.create_button') }}
        </UButton>
      </div>
    </template>
  </UModal>
</template>
