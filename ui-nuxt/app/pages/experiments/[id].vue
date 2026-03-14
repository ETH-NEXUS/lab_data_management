<script setup lang="ts">
import { computed, ref } from 'vue'
import { useQueryClient } from '@tanstack/vue-query'
import { useExperimentQuery } from '~/composables/useExperimentQuery'
import { useProjectQuery } from '~/composables/useProjectsQuery'
import { EXPERIMENTS_QUERY_KEY, getExperimentQueryKey } from '~/types/experiments'
import { getProjectQueryKey, PROJECTS_QUERY_KEY } from '~/types/projects'
import { formatDateTime } from '~/utils/dateTime'
import { getErrorMessage } from '~/utils/errors'
import ExperimentBarcodeSpecificationModal from '~/components/experiments/ExperimentBarcodeSpecificationModal.vue'
import ExperimentBarcodeSpecificationsCard from '~/components/experiments/ExperimentBarcodeSpecificationsCard.vue'
import ExperimentHeader from '~/components/experiments/ExperimentHeader.vue'
import ExperimentInfoCard from '~/components/experiments/ExperimentInfoCard.vue'
import ExperimentEditModal from '~/components/experiments/ExperimentEditModal.vue'

const route = useRoute()
const { t } = useI18n()
const queryClient = useQueryClient()

const experimentId = computed(() => Number(route.params.id))
const hasValidExperimentId = computed(() => Number.isInteger(experimentId.value) && experimentId.value > 0)

const experimentQuery = useExperimentQuery(experimentId)
const experiment = computed(() => experimentQuery.data.value ?? null)

const projectId = computed(() => experiment.value?.project ?? 0)
const projectQuery = useProjectQuery(projectId)
const projectName = computed(() => projectQuery.data.value?.name ?? null)

const isEditModalOpen = ref(false)
const isAddBarcodeSpecificationModalOpen = ref(false)
const editField = ref<'name' | 'description'>('name')
const editInitialValue = ref('')

const experimentErrorMessage = computed(() => {
  const err = experimentQuery.error.value
  if (!err) return null
  return getErrorMessage(err)
})

/**
 * Invalidates experiment and project-related query caches after experiment updates.
 *
 * Invalidated key examples:
 * - `['experiment', '21']` for the open experiment detail page
 * - `['experiments']` for global experiment list views
 * - `['projects']` for navigation tree where experiments are nested under projects
 * - `['project', '5']` for a loaded parent project detail page
 */
const invalidateExperimentQueries = async () => {
  if (!hasValidExperimentId.value) return

  const invalidations: Array<Promise<void>> = [
    queryClient.invalidateQueries({ queryKey: getExperimentQueryKey(experimentId.value) }),
    queryClient.invalidateQueries({ queryKey: EXPERIMENTS_QUERY_KEY }),
    queryClient.invalidateQueries({ queryKey: PROJECTS_QUERY_KEY }),
  ]

  if (projectId.value > 0) {
    invalidations.push(queryClient.invalidateQueries({ queryKey: getProjectQueryKey(projectId.value) }))
  }

  await Promise.all(invalidations)
}

const openEditModal = (field: 'name' | 'description') => {
  if (!experiment.value) return

  editField.value = field
  editInitialValue.value = field === 'name' ? experiment.value.name : (experiment.value.description ?? '')
  isEditModalOpen.value = true
}

const onExperimentSaved = async () => {
  await invalidateExperimentQueries()
}

const openAddBarcodeSpecificationModal = () => {
  isAddBarcodeSpecificationModalOpen.value = true
}

const onBarcodeSpecificationCreated = async () => {
  await invalidateExperimentQueries()
}

const onBarcodeSpecificationDeleted = async () => {
  await invalidateExperimentQueries()
}
</script>

<template>
  <section class="space-y-5 p-6">
    <p v-if="experimentQuery.isPending.value" class="text-sm text-slate-500">{{ t('experiments.page.loading') }}</p>
    <p v-else-if="experimentErrorMessage" class="rounded-md border border-red-200 bg-red-50 p-3 text-sm text-red-700">
      {{ experimentErrorMessage }}
    </p>
    <p v-else-if="!experiment" class="text-sm text-slate-500">{{ t('experiments.page.not_found') }}</p>

    <template v-else>
      <ExperimentHeader :experiment="experiment" :project-name="projectName" @edit-name="openEditModal('name')" />

      <ExperimentInfoCard
        :experiment="experiment"
        :created-at-label="formatDateTime(experiment.created_at)"
        @edit-description="openEditModal('description')"
      />

      <ExperimentBarcodeSpecificationsCard
        :experiment="experiment"
        @add-specification="openAddBarcodeSpecificationModal"
        @deleted-specification="onBarcodeSpecificationDeleted"
      />

      <ExperimentEditModal
        v-model:open="isEditModalOpen"
        :experiment-id="experiment.id"
        :field="editField"
        :initial-value="editInitialValue"
        @saved="onExperimentSaved"
      />

      <ExperimentBarcodeSpecificationModal
        v-model:open="isAddBarcodeSpecificationModalOpen"
        :experiment-id="experiment.id"
        @created="onBarcodeSpecificationCreated"
      />
    </template>
  </section>
</template>
