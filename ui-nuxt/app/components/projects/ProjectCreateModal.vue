<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import { useHarvestStore } from '~/stores/harvest'
import { useProjectStore } from '~/stores/projects'
import type { CreateProjectPayload, Project } from '~/types/projects'
import { getErrorMessage } from '~/utils/errors'

const props = defineProps<{
  open: boolean
}>()

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
  (e: 'created', project: Project): void
}>()

const toast = useToast()
const { t } = useI18n()
const harvestStore = useHarvestStore()
const projectStore = useProjectStore()

const useCustomName = ref(false)
const customProjectName = ref('')
const selectedHarvestProjectId = ref<number | null>(null)

/**
 * Loads Harvest projects when the modal opens and resets local modal state.
 *
 * Local state example after opening:
 * - `{ useCustomName: false, customProjectName: '', selectedHarvestProjectId: null }`
 */
watch(
  () => props.open,
  async (isOpen) => {
    if (!isOpen) return

    useCustomName.value = false
    customProjectName.value = ''
    selectedHarvestProjectId.value = null

    try {
      await harvestStore.initialize()
    } catch {
      // Keep modal usable with manual name even when Harvest loading fails.
    }
  },
  { immediate: true },
)

const hasHarvestProjects = computed(() => harvestStore.harvestProjects.length > 0)

const selectedHarvestProject = computed(
  () => harvestStore.harvestProjects.find((project) => project.id === selectedHarvestProjectId.value) ?? null,
)

/**
 * Decides which project name should be sent to the backend.
 *
 * Output examples:
 * - `'Screening 2026'` (picked from Harvest)
 * - `'My custom project'` (typed manually)
 * - `''` (invalid, save disabled)
 */
const resolvedProjectName = computed(() => {
  if (useCustomName.value || !hasHarvestProjects.value) {
    return customProjectName.value.trim()
  }

  return selectedHarvestProject.value?.name?.trim() ?? ''
})

const canSave = computed(() => resolvedProjectName.value.length > 0 && !projectStore.isCreatingProject)

const close = () => emit('update:open', false)

/**
 * Creates payload that mirrors the old behavior:
 * - custom name => no Harvest relation
 * - Harvest selection => set harvest_id + harvest_notes
 *
 * Payload examples:
 * - `{ name: 'Manual', harvest_id: null, harvest_notes: null }`
 * - `{ name: 'H-Project', harvest_id: 7, harvest_notes: 'Harvest notes' }`
 */
const buildCreatePayload = (): CreateProjectPayload => {
  if (useCustomName.value || !hasHarvestProjects.value) {
    return {
      name: resolvedProjectName.value,
      harvest_id: null,
      harvest_notes: null,
    }
  }

  return {
    name: resolvedProjectName.value,
    harvest_id: selectedHarvestProject.value?.id ?? null,
    harvest_notes: selectedHarvestProject.value?.notes ?? null,
  }
}

const refreshHarvestProjects = async () => {
  try {
    await harvestStore.getHarvestProjects()
  } catch (err: unknown) {
    toast.add({
      title: t('projects.create_modal.harvest_refresh_failed'),
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  }
}

const save = async () => {
  if (!canSave.value) return

  try {
    const payload = buildCreatePayload()
    const project = await projectStore.createProject(payload)

    emit('created', project)
    close()

    toast.add({
      title: t('projects.create_modal.created'),
      description: project.name,
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: t('projects.create_modal.create_failed'),
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
    :title="t('projects.create_modal.title')"
    :description="t('projects.create_modal.description')"
    :dismissible="!projectStore.isCreatingProject"
    class="w-full sm:max-w-3xl"
    @update:open="emit('update:open', $event)"
  >
    <template #body>
      <div class="space-y-4">
        <div v-if="!useCustomName && hasHarvestProjects" class="space-y-2">
          <p class="text-sm font-semibold text-slate-700">{{ t('projects.create_modal.source_harvest_label') }}</p>

          <USelectMenu
            v-model="selectedHarvestProjectId"
            :items="harvestStore.harvestProjects"
            value-key="id"
            label-key="name"
            :placeholder="t('projects.create_modal.harvest_placeholder')"
            class="w-full"
          />

          <UButton
            size="xs"
            variant="soft"
            color="secondary"
            icon="i-heroicons-arrow-path"
            :loading="harvestStore.isLoading"
            @click="refreshHarvestProjects"
          >
            {{ t('projects.create_modal.refresh_harvest') }}
          </UButton>
        </div>

        <div v-if="useCustomName || !hasHarvestProjects" class="space-y-2">
          <p class="text-sm font-semibold text-slate-700">{{ t('projects.create_modal.custom_name_label') }}</p>
          <UInput
            v-model="customProjectName"
            class="w-full"
            :placeholder="t('projects.create_modal.custom_name_placeholder')"
            autofocus
          />
        </div>

        <div v-if="hasHarvestProjects">
          <UButton
            size="xs"
            variant="ghost"
            color="neutral"
            :icon="useCustomName ? 'i-heroicons-list-bullet' : 'i-heroicons-pencil-square'"
            @click="useCustomName = !useCustomName"
          >
            {{
              useCustomName
                ? t('projects.create_modal.toggle_use_harvest')
                : t('projects.create_modal.toggle_use_custom')
            }}
          </UButton>
        </div>
      </div>
    </template>

    <template #footer>
      <div class="flex w-full justify-end gap-2">
        <UButton color="neutral" variant="soft" :disabled="projectStore.isCreatingProject" @click="close">
          {{ t('common.actions.cancel') }}
        </UButton>

        <UButton
          color="primary"
          variant="solid"
          :loading="projectStore.isCreatingProject"
          :disabled="!canSave"
          @click="save"
        >
          {{ t('projects.create_modal.create_button') }}
        </UButton>
      </div>
    </template>
  </UModal>
</template>
