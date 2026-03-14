<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import BaseField from '~/components/common/BaseField.vue'
import WavesModalWrapper from '~/components/common/WavesModalWrapper.vue'
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

const textAreaInputClassName =
  'w-full px-4 py-3 bg-white outline-none ring-offset-0 focus:ring-2 focus:ring-lime-500 shadow rounded-2xl'

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
  <WavesModalWrapper
    :open="props.open"
    :title="title"
    :description="description"
    :dismissible="!projectStore.isUpdatingProject"
    modal-class="w-full sm:max-w-3xl"
    body-container-class="w-full max-w-2xl px-8 pt-10"
    @update:open="emit('update:open', $event)"
  >
    <template #body>
      <div class="space-y-6">
        <BaseField
          v-if="props.field === 'name'"
          v-model="localValue"
          :label="t('projects.edit_modal.field_name')"
          :placeholder="t('projects.edit_modal.placeholder_name')"
          :autofocus="true"
        />
        <BaseField
          v-else
          v-model="localValue"
          :label="t('projects.edit_modal.field_description')"
          :placeholder="t('projects.edit_modal.placeholder_description')"
          :autofocus="true"
          :multiline="true"
          :rows="5"
          :input-class="textAreaInputClassName"
        />
      </div>
    </template>

    <template #footer>
      <BaseButton
        :label="t('common.actions.cancel')"
        :on-click="close"
        variant="secondary"
        size="sm"
        width="auto"
        :disabled="projectStore.isUpdatingProject"
      />
      <BaseButton
        :label="t('common.actions.save')"
        :on-click="save"
        variant="primary"
        size="sm"
        width="auto"
        :loading="projectStore.isUpdatingProject"
        :disabled="!canSave"
      />
    </template>
  </WavesModalWrapper>
</template>
