<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import { useProjectStore } from '~/stores/projects'
import type { ProjectPayload } from '~/types/projects'
import { getErrorMessage } from '~/utils/errors'

const props = defineProps<{
  open: boolean
  projectId: number
  field: 'name' | 'description'
  initialValue: string
}>()

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
  (e: 'saved'): void
}>()

const toast = useToast()
const { t } = useI18n()
const projectStore = useProjectStore()
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
  return props.field === 'name' ? t('projects.edit_modal.title_name') : t('projects.edit_modal.title_description')
})

const description = computed(() => {
  return props.field === 'name'
    ? t('projects.edit_modal.description_name')
    : t('projects.edit_modal.description_description')
})

const canSave = computed(() => {
  if (props.field === 'name') {
    return localValue.value.trim().length > 0
  }

  return true
})

const close = () => emit('update:open', false)

const save = async () => {
  if (!canSave.value) return

  const payload = { [props.field]: localValue.value } as ProjectPayload

  try {
    await projectStore.updateProject(props.projectId, payload)
    emit('saved')
    close()
    toast.add({
      title: t('projects.edit_modal.updated'),
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: t('projects.edit_modal.update_failed'),
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
    :dismissible="!projectStore.isUpdatingProject"
    class="w-full sm:max-w-3xl"
    @update:open="emit('update:open', $event)"
  >
    <template #body>
      <div class="space-y-2">
        <p class="text-sm text-slate-600">
          {{
            props.field === 'name' ? t('projects.edit_modal.field_name') : t('projects.edit_modal.field_description')
          }}
        </p>

        <UInput
          v-if="props.field === 'name'"
          v-model="localValue"
          class="w-full"
          :placeholder="t('projects.edit_modal.placeholder_name')"
          autofocus
        />
        <UTextarea
          v-else
          v-model="localValue"
          class="w-full"
          :placeholder="t('projects.edit_modal.placeholder_description')"
          :rows="5"
          autofocus
        />
      </div>
    </template>

    <template #footer>
      <div class="flex w-full justify-end gap-2">
        <UButton color="neutral" variant="soft" :disabled="projectStore.isUpdatingProject" @click="close">
          {{ t('common.actions.cancel') }}
        </UButton>
        <UButton
          color="primary"
          variant="solid"
          :loading="projectStore.isUpdatingProject"
          :disabled="!canSave"
          @click="save"
        >
          {{ t('common.actions.save') }}
        </UButton>
      </div>
    </template>
  </UModal>
</template>
