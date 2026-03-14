<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import BaseField from '~/components/common/BaseField.vue'
import WavesModalWrapper from '~/components/common/WavesModalWrapper.vue'
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
 * Keeps native select value in sync with numeric project id state.
 *
 * Data examples:
 * - state `{ selectedHarvestProjectId: 8 }` => select model `'8'`
 * - select model `''` => state `{ selectedHarvestProjectId: null }`
 */
const selectedHarvestProjectIdText = computed({
  get: (): string => {
    if (selectedHarvestProjectId.value === null) return ''
    return String(selectedHarvestProjectId.value)
  },
  set: (value: string) => {
    if (value.trim() === '') {
      selectedHarvestProjectId.value = null
      return
    }

    const parsedValue = Number(value)
    selectedHarvestProjectId.value = Number.isNaN(parsedValue) ? null : parsedValue
  },
})

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

const customNameInputClassName =
  'w-full rounded-full border border-slate-300 bg-white px-4 py-3 outline-none ring-offset-0 focus:ring-2 focus:ring-lime-500'

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
  <WavesModalWrapper
    :open="props.open"
    :title="t('projects.create_modal.title')"
    :description="t('projects.create_modal.description')"
    :dismissible="!projectStore.isCreatingProject"
    modal-class="w-full sm:max-w-4xl"
    body-container-class="w-full max-w-2xl px-8 pt-10"
    @update:open="emit('update:open', $event)"
  >
    <template #body>
      <div class="space-y-3">
        <div v-if="!useCustomName && hasHarvestProjects" class="space-y-2">
          <label class="mb-1 block pl-4 text-sm font-medium">
            {{ t('projects.create_modal.source_harvest_label') }}
          </label>
          <div class="relative w-full rounded-full bg-white">
            <span class="pointer-events-none absolute top-1/2 right-0 mr-5 -translate-y-1/2 text-slate-500">
              <svg width="24" height="24" viewBox="0 0 20 20" fill="none" xmlns="http://www.w3.org/2000/svg">
                <path d="M6.4 8.2L10 11.8L13.6 8.2" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" />
              </svg>
            </span>
            <select
              v-model="selectedHarvestProjectIdText"
              class="relative w-full cursor-pointer appearance-none rounded-full border border-slate-300 bg-transparent px-4 py-3 pr-14 text-slate-600 ring-offset-0 transition outline-none focus:ring-2 focus:ring-lime-500"
            >
              <option value="">{{ t('projects.create_modal.harvest_placeholder') }}</option>
              <option v-for="project in harvestStore.harvestProjects" :key="project.id" :value="String(project.id)">
                {{ project.name }}
              </option>
            </select>
          </div>

          <button
            type="button"
            class="mt-1 inline-block cursor-pointer text-sm font-medium text-teal-800 transition hover:text-teal-600 hover:underline disabled:cursor-not-allowed disabled:opacity-60"
            :disabled="harvestStore.isLoading"
            @click="refreshHarvestProjects"
          >
            {{ t('projects.create_modal.refresh_harvest') }}
          </button>
        </div>

        <div v-if="useCustomName || !hasHarvestProjects" class="space-y-2">
          <BaseField
            v-model="customProjectName"
            field-class="mb-6"
            :label="t('projects.create_modal.custom_name_label')"
            :placeholder="t('projects.create_modal.custom_name_placeholder')"
            :input-class="customNameInputClassName"
            :autofocus="true"
          />
        </div>

        <div v-if="hasHarvestProjects">
          <BaseButton
            :label="
              useCustomName
                ? t('projects.create_modal.toggle_use_harvest')
                : t('projects.create_modal.toggle_use_custom')
            "
            :on-click="
              () => {
                useCustomName = !useCustomName
              }
            "
            variant="accent"
            size="sm"
            width="auto"
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
        :disabled="projectStore.isCreatingProject"
      />

      <BaseButton
        :label="t('projects.create_modal.create_button')"
        :on-click="save"
        variant="primary"
        size="sm"
        width="auto"
        :loading="projectStore.isCreatingProject"
        :disabled="!canSave"
      />
    </template>
  </WavesModalWrapper>
</template>
