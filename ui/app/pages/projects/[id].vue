<script setup lang="ts">
import { computed, ref } from 'vue'
import { useQueryClient } from '@tanstack/vue-query'
import { useProjectQuery } from '~/composables/useProjectsQuery'
import { useHarvestStore } from '~/stores/harvest'
import { getProjectQueryKey, PROJECTS_QUERY_KEY } from '~/types/projects'
import { EXPERIMENTS_QUERY_KEY } from '~/types/experiments'
import { formatDateTime } from '~/utils/dateTime'
import { getErrorMessage } from '~/utils/errors'
import ProjectHeader from 'components/projects/ProjectHeader.vue'
import ProjectOverviewCard from 'components/projects/ProjectOverviewCard.vue'
import ProjectExperimentsCard from 'components/projects/ProjectExperimentsCard.vue'
import ProjectEditModal from 'components/projects/ProjectEditModal.vue'

const route = useRoute()
const toast = useToast()
const { t } = useI18n()
const queryClient = useQueryClient()
const harvestStore = useHarvestStore()

const projectId = computed(() => Number(route.params.id))
const hasValidProjectId = computed(() => Number.isInteger(projectId.value) && projectId.value > 0)

const projectQuery = useProjectQuery(projectId)

const project = computed(() => projectQuery.data.value ?? null)
const experimentsCount = computed(() => project.value?.experiments?.length ?? 0)

const isEditModalOpen = ref(false)
const editField = ref<'name' | 'description'>('name')
const editInitialValue = ref('')

const projectErrorMessage = computed(() => {
  const err = projectQuery.error.value
  if (!err) return null
  return getErrorMessage(err)
})

/**
 * Marks both related caches as stale after a project mutation so Vue Query refetches fresh data.
 *
 * Invalidated query examples:
 * - `{ queryKey: ['project', '12'] }` for the current project detail page
 * - `{ queryKey: ['projects'] }` for the navigation/list that shows all projects
 */
const invalidateProjectQueries = async () => {
  if (!hasValidProjectId.value) return

  await Promise.all([
    queryClient.invalidateQueries({ queryKey: getProjectQueryKey(projectId.value) }),
    queryClient.invalidateQueries({ queryKey: PROJECTS_QUERY_KEY }),
    queryClient.invalidateQueries({ queryKey: EXPERIMENTS_QUERY_KEY }),
  ])
}

const openEditModal = (field: 'name' | 'description') => {
  if (!project.value) return

  editField.value = field
  editInitialValue.value = field === 'name' ? project.value.name : (project.value.description ?? '')
  isEditModalOpen.value = true
}

const onProjectSaved = async () => {
  await invalidateProjectQueries()
}

const onExperimentsMerged = async () => {
  await invalidateProjectQueries()
}

const updateHarvestInfo = async () => {
  if (!project.value) return

  try {
    const resp = await harvestStore.updateHarvestInfo(project.value.id)
    await invalidateProjectQueries()

    toast.add({
      title: resp.success ? t('projects.page.harvest_update_success') : t('projects.page.harvest_update_requested'),
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: t('projects.page.harvest_update_failed'),
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  }
}
</script>

<template>
  <section class="space-y-5 p-6">
    <p v-if="projectQuery.isPending.value" class="text-sm text-slate-500">{{ t('projects.page.loading') }}</p>
    <p v-else-if="projectErrorMessage" class="rounded-md border border-red-200 bg-red-50 p-3 text-sm text-red-700">
      {{ projectErrorMessage }}
    </p>
    <p v-else-if="!project" class="text-sm text-slate-500">{{ t('projects.page.not_found') }}</p>

    <template v-else>
      <ProjectHeader :project="project" @edit-name="openEditModal('name')" />

      <ProjectOverviewCard
        :project="project"
        :experiments-count="experimentsCount"
        :created-at-label="formatDateTime(project.created_at)"
        :harvest-loading="harvestStore.isUpdatingInfo"
        @edit-description="openEditModal('description')"
        @update-harvest="updateHarvestInfo"
      />

      <ProjectExperimentsCard :project="project" @moved="onExperimentsMerged" />

      <ProjectEditModal
        v-model:open="isEditModalOpen"
        :project-id="project.id"
        :field="editField"
        :initial-value="editInitialValue"
        @saved="onProjectSaved"
      />
    </template>
  </section>
</template>
