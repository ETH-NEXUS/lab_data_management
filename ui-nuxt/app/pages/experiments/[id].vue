<script setup lang="ts">
import { computed, ref, watch } from 'vue'
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
import ExperimentHeatmap from '~/components/experiments/ExperimentHeatmap.vue'
import ExperimentMeasurementCalculatorModal from '~/components/experiments/ExperimentMeasurementCalculatorModal.vue'
import ExperimentGenerateReportModal from '~/components/experiments/ExperimentGenerateReportModal.vue'
import { usePlateViewStore } from '~/stores/plateView'
import { useReportStore } from '~/stores/reports'

const route = useRoute()
const { t } = useI18n()
const toast = useToast()
const queryClient = useQueryClient()
const reportStore = useReportStore()
const plateViewStore = usePlateViewStore()

const experimentId = computed(() => Number(route.params.id))
const hasValidExperimentId = computed(() => Number.isInteger(experimentId.value) && experimentId.value > 0)

const experimentQuery = useExperimentQuery(experimentId)
const experiment = computed(() => experimentQuery.data.value ?? null)

const projectId = computed(() => experiment.value?.project ?? 0)
const projectQuery = useProjectQuery(projectId)
const projectName = computed(() => projectQuery.data.value?.name ?? null)
const canAddExperimentData = computed(() => {
  if (!experiment.value) return false
  if (!experiment.value.details) return false
  if (!experiment.value.details.measurement_labels) return false
  return experiment.value.details.measurement_labels.length > 0
})
const canAddMeasurement = computed(() => {
  if (!experiment.value) return false
  if (!experiment.value.available_measurement_labels) return false
  return experiment.value.available_measurement_labels.length > 0
})
const canShowResults = computed(() => {
  if (!experiment.value) return false
  if (!experiment.value.details) return false
  if (!experiment.value.details.measurement_labels) return false
  return experiment.value.details.measurement_labels.length > 0
})
const canGenerateReport = computed(() => {
  if (!experiment.value) return false
  if (!experiment.value.available_measurement_labels) return false
  return experiment.value.available_measurement_labels.length > 0
})

const isEditModalOpen = ref(false)
const isAddBarcodeSpecificationModalOpen = ref(false)
const isMeasurementCalculatorModalOpen = ref(false)
const isGenerateReportModalOpen = ref(false)
const isResultsExpanded = ref(false)
const editField = ref<'name' | 'description'>('name')
const editInitialValue = ref('')
const selectedPosControl = ref<string | null>(null)
const selectedNegControl = ref<string | null>(null)

const experimentErrorMessage = computed(() => {
  const err = experimentQuery.error.value
  if (!err) return null
  return getErrorMessage(err)
})

watch(
  () => experiment.value?.name,
  async (experimentName) => {
    if (!experimentName) {
      reportStore.outputNotebooks = []
      return
    }

    await reportStore.fetchOutputNotebooks(experimentName)
  },
  { immediate: true },
)

watch(
  () => experiment.value?.id,
  (currentExperimentId) => {
    if (!currentExperimentId) {
      selectedPosControl.value = null
      selectedNegControl.value = null
      return
    }

    const storedControls = plateViewStore.getExperimentControls(currentExperimentId)
    if (!storedControls) {
      selectedPosControl.value = null
      selectedNegControl.value = null
      return
    }

    selectedPosControl.value = storedControls.pos
    selectedNegControl.value = storedControls.neg
  },
  { immediate: true },
)

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

const openAddExperimentDataPage = async () => {
  if (!experiment.value) return
  await navigateTo(`/add_data/${experiment.value.id}`)
}

const openMeasurementCalculatorModal = () => {
  isMeasurementCalculatorModalOpen.value = true
}

const onMeasurementCalculated = async () => {
  await invalidateExperimentQueries()
}

const openGenerateReportModal = () => {
  isGenerateReportModalOpen.value = true
}

const onReportGenerated = async () => {
  if (!experiment.value) {
    return
  }

  await reportStore.fetchOutputNotebooks(experiment.value.name)
}

/**
 * Shows short filename from one full report path.
 *
 * Accepted input example:
 * - `'/notebooks/output/dose_response_report.pdf'`
 *
 * Returned data example:
 * - `'dose_response_report.pdf'`
 */
const getReportLabel = (reportPath: string): string => {
  return reportStore.getFileNameFromPath(reportPath)
}

const downloadReport = async (reportPath: string) => {
  try {
    await reportStore.downloadReport(reportPath)
  } catch (err: unknown) {
    toast.add({
      title: t('experiments.reports.download_failed'),
      description: getErrorMessage(err),
      color: 'error',
      duration: 4500,
    })
  }
}

const onResultsToggle = (event: Event) => {
  const target = event.target as HTMLDetailsElement | null
  if (!target) return
  isResultsExpanded.value = target.open
}

const handleUpdateControls = (data: { pos: string | null; neg: string | null }) => {
  selectedPosControl.value = data.pos
  selectedNegControl.value = data.neg
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

      <section v-if="canAddExperimentData || canAddMeasurement || canGenerateReport" class="mx-auto w-[80%]">
        <div class="flex flex-wrap gap-2">
          <UButton
            v-if="canAddExperimentData"
            color="secondary"
            variant="solid"
            icon="i-heroicons-table-cells"
            :label="t('experiments.page.add_experiment_data')"
            @click="openAddExperimentDataPage"
          />
          <UButton
            v-if="canAddMeasurement"
            color="secondary"
            variant="outline"
            icon="i-heroicons-calculator"
            :label="t('experiments.page.add_measurement')"
            @click="openMeasurementCalculatorModal"
          />
          <UButton
            v-if="canGenerateReport"
            color="secondary"
            variant="outline"
            icon="i-heroicons-document-arrow-down"
            :label="t('experiments.page.generate_report')"
            @click="openGenerateReportModal"
          />
        </div>
      </section>

      <ExperimentBarcodeSpecificationsCard
        :experiment="experiment"
        @add-specification="openAddBarcodeSpecificationModal"
        @deleted-specification="onBarcodeSpecificationDeleted"
      />

      <UCard
        v-if="reportStore.isLoadingOutputNotebooks || reportStore.outputNotebooks.length > 0"
        class="mx-auto w-[80%]"
        :ui="{
          root: 'core-card divide-y divide-slate-200/70',
        }"
      >
        <template #header>
          <p class="font-semibold">{{ t('experiments.reports.available_reports') }}</p>
        </template>

        <p v-if="reportStore.isLoadingOutputNotebooks" class="text-sm text-slate-600">
          {{ t('experiments.reports.loading_reports') }}
        </p>

        <div v-else class="space-y-2">
          <button
            v-for="reportPath in reportStore.outputNotebooks"
            :key="`report-${reportPath}`"
            type="button"
            class="cursor-pointer text-sm text-blue-700 hover:text-blue-800 hover:underline"
            @click="downloadReport(reportPath)"
          >
            {{ getReportLabel(reportPath) }}
          </button>
        </div>
      </UCard>

      <UCard
        v-if="canShowResults"
        class="mx-auto w-[80%]"
        :ui="{
          root: 'core-card divide-y divide-slate-200/70',
        }"
      >
        <details class="group rounded-xl border border-slate-200 bg-white/80 p-3" @toggle="onResultsToggle">
          <summary class="flex cursor-pointer list-none items-center justify-between gap-3">
            <span class="text-sm font-semibold text-slate-800">{{ t('experiments.results.show_results') }}</span>
            <UIcon
              name="i-heroicons-chevron-right"
              class="size-5 shrink-0 text-slate-500 transition-transform duration-200 group-open:rotate-90"
            />
          </summary>

          <div v-if="isResultsExpanded" class="mt-4">
            <ExperimentHeatmap
              :experiment-id="experiment.id"
              :timestamps="experiment.details.measurement_timestamps"
              :available-measurement-labels="experiment.details.measurement_labels"
              :overall-stats="experiment.details.overall_stats"
              @update-controls="handleUpdateControls"
            />
          </div>
        </details>
      </UCard>

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

      <ExperimentMeasurementCalculatorModal
        v-model:open="isMeasurementCalculatorModalOpen"
        :experiment-id="experiment.id"
        :labels="experiment.available_measurement_labels"
        :measurement-timestamps="experiment.details.measurement_timestamps"
        @calculated="onMeasurementCalculated"
      />

      <ExperimentGenerateReportModal
        v-model:open="isGenerateReportModalOpen"
        :experiment-name="experiment.name"
        :labels="experiment.available_measurement_labels"
        :selected-pos="selectedPosControl"
        :selected-neg="selectedNegControl"
        @generated="onReportGenerated"
      />
    </template>
  </section>
</template>
