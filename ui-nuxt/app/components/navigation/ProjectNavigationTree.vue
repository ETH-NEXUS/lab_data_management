<script setup lang="ts">
import { computed, ref } from 'vue'
import { useQueryClient } from '@tanstack/vue-query'
import { useProjectsQuery } from '~/composables/useProjectsQuery'
import type { NavigationTreeNode } from '~/types/navigation'
import { EXPERIMENTS_QUERY_KEY, type Experiment } from '~/types/experiments'
import { PROJECTS_QUERY_KEY, type Project } from '~/types/projects'
import ExperimentCreateModal from '~/components/experiments/ExperimentCreateModal.vue'
import ProjectCreateModal from '~/components/projects/ProjectCreateModal.vue'

const props = withDefaults(defineProps<{ filter?: string }>(), {
  filter: '',
})

const { t } = useI18n()
const queryClient = useQueryClient()
const isCreateProjectModalOpen = ref(false)
const isCreateExperimentModalOpen = ref(false)
const selectedProjectIdForNewExperiment = ref<number | null>(null)
const selectedProjectNameForNewExperiment = ref<string | null>(null)
const projectsQuery = useProjectsQuery()

/**
 * Context menu for the "Projects" root node.
 *
 * Returned menu shape example:
 * - `[[{ label: 'Create project', icon: 'i-heroicons-plus', onSelect: () => {...} }]]`
 */
const projectsRootContextMenuItems = computed(() => [
  [
    {
      label: t('navigation.projects.actions.create_project'),
      icon: 'i-heroicons-plus-circle',
      onSelect: () => {
        isCreateProjectModalOpen.value = true
      },
    },
  ],
])

/**
 * Maps one project's experiments into child tree nodes.
 *
 * Input example:
 * - `{ id: 5, name: 'Screening', experiments: [{ id: 21, name: 'Dose response', plates: [{ id: 7, barcode: 'A001', dimension: { id: 1, name: '384', cols: 24, rows: 16 } }] }] }`
 *
 * Returned data example:
 * - `[{ label: 'Dose response', icon: 'i-heroicons-beaker', children: [{ label: 'A001 (384)', icon: 'i-heroicons-squares-2x2' }] }]`
 */
const formatPlateDimensionLabel = (dimension: unknown): string => {
  if (typeof dimension === 'string' || typeof dimension === 'number') {
    return String(dimension)
  }

  if (typeof dimension === 'object' && dimension !== null && 'name' in dimension) {
    const name = (dimension as { name?: unknown }).name
    if (typeof name === 'string' && name.trim() !== '') {
      return name
    }
  }

  return t('navigation.projects.no_dimension')
}

/**
 * Maps one experiment's plates into third-level tree nodes.
 *
 * Input example:
 * - `{ id: 21, name: 'Dose response', plates: [{ id: 11, barcode: 'B001', dimension: { id: 2, name: '96', cols: 12, rows: 8 } }, { id: 10, barcode: 'A001', dimension: { id: 1, name: '384', cols: 24, rows: 16 } }] }`
 *
 * Returned data example:
 * - `[{ label: 'A001 (384)', icon: 'i-heroicons-squares-2x2' }, { label: 'B001 (96)', icon: 'i-heroicons-squares-2x2' }]`
 */
const mapPlateNodes = (experiment: Experiment): NavigationTreeNode[] => {
  const sortedPlates = [...(experiment.plates ?? [])].sort((a, b) => a.barcode.localeCompare(b.barcode))

  return sortedPlates.map((plate) => ({
    label: `${plate.barcode || t('navigation.projects.plate_fallback', { id: plate.id })} (${formatPlateDimensionLabel(plate.dimension)})`,
    icon: 'i-heroicons-squares-2x2',
    /**
     * Navigates to the migrated plate page.
     *
     * Route param examples:
     * - barcode present: `{ id: 'ABC-0001' }` -> `/plates/ABC-0001`
     * - barcode missing: `{ id: 17 }` -> `/plates/17`
     */
    onSelect: () => navigateTo(`/plates/${encodeURIComponent(plate.barcode || String(plate.id))}`),
  }))
}

const mapExperimentNodes = (project: Project): NavigationTreeNode[] => {
  return (project.experiments ?? []).map((experiment) => ({
    label: experiment.name || t('navigation.projects.experiment_fallback', { id: experiment.id }),
    icon: 'i-heroicons-beaker',
    onSelect: () => navigateTo(`/experiments/${experiment.id}`),
    children: mapPlateNodes(experiment),
  }))
}

const projectNodes = computed<NavigationTreeNode[]>(() =>
  (projectsQuery.data.value ?? []).map((project) => ({
    label: project.name || t('navigation.projects.project_fallback', { id: project.id }),
    slot: 'project-node',
    projectId: project.id,
    projectName: project.name || t('navigation.projects.project_fallback', { id: project.id }),
    icon: 'i-heroicons-document-text',
    onSelect: () => navigateTo(`/projects/${project.id}`),
    children: mapExperimentNodes(project),
  })),
)

const rawItems = computed<NavigationTreeNode[]>(() => [
  {
    label: t('navigation.projects.root_label'),
    slot: 'projects-root',
    icon: 'i-heroicons-folder-open',
    defaultExpanded: true,
    children: projectNodes.value,
  },
])

const filterTree = (items: NavigationTreeNode[], query: string): NavigationTreeNode[] => {
  const q = query.trim().toLowerCase()
  if (!q) return items

  return items
    .map((item) => {
      const children = item.children ? filterTree(item.children, q) : []
      const selfMatch = item.label.toLowerCase().includes(q)

      if (!selfMatch && children.length === 0) return null

      return {
        ...item,
        children,
        defaultExpanded: true,
      }
    })
    .filter((item): item is NavigationTreeNode => item !== null)
}

const items = computed<NavigationTreeNode[]>(() => filterTree(rawItems.value, props.filter))

const projectsErrorMessage = computed(() => {
  const err = projectsQuery.error.value
  if (!err) return null
  return err instanceof Error ? err.message : t('navigation.projects.error')
})

/**
 * Opens the create experiment modal for one project.
 *
 * Input example:
 * - `{ projectId: 5, projectName: 'Screening' }`
 */
const openCreateExperimentModal = (projectId: number, projectName: string) => {
  selectedProjectIdForNewExperiment.value = projectId
  selectedProjectNameForNewExperiment.value = projectName
  isCreateExperimentModalOpen.value = true
}

/**
 * Context menu entries for project node actions.
 *
 * Returned menu shape example:
 * - `[[{ label: 'Create experiment', onSelect: () => {...} }]]`
 */
const getProjectContextMenuItems = (projectId: number, projectName: string) => [
  [
    {
      label: t('navigation.projects.actions.create_experiment'),
      icon: 'i-heroicons-plus-circle',
      onSelect: () => openCreateExperimentModal(projectId, projectName),
    },
  ],
]

/**
 * Runtime-safe extraction of `label` from UTree slot item payload.
 *
 * Accepted payload examples:
 * - `{ label: 'Projects' }`
 * - `{ label: 'Screening 2026', projectId: 5, projectName: 'Screening 2026' }`
 *
 * Returned value examples:
 * - `'Projects'`
 * - `'Screening 2026'`
 * - `''` (when label is missing)
 */
const getSlotItemLabel = (item: unknown): string => {
  if (typeof item !== 'object' || item === null || !('label' in item)) return ''

  const label = (item as { label?: unknown }).label
  return typeof label === 'string' ? label : ''
}

/**
 * Builds project-node context menu items from UTree slot item payload.
 *
 * Accepted payload example:
 * - `{ label: 'Screening 2026', projectId: 5, projectName: 'Screening 2026' }`
 *
 * Returned menu shape example:
 * - `[[{ label: 'Create experiment', icon: 'i-heroicons-plus-circle', onSelect: () => void }]]`
 *
 * Fallback behavior:
 * - returns `[]` when payload does not contain a valid project id/name.
 */
const getProjectContextMenuItemsFromSlotItem = (item: unknown) => {
  if (typeof item !== 'object' || item === null) return []
  if (!('projectId' in item) || !('projectName' in item)) return []

  const projectId = (item as { projectId?: unknown }).projectId
  const projectName = (item as { projectName?: unknown }).projectName

  if (typeof projectId !== 'number' || typeof projectName !== 'string') return []
  return getProjectContextMenuItems(projectId, projectName)
}

/**
 * Refreshes projects list cache after creating a new project.
 *
 * Input example:
 * - `{ id: 44, name: 'Screening 2026' }`
 */
const onProjectCreated = async (_project: Project) => {
  await queryClient.invalidateQueries({ queryKey: PROJECTS_QUERY_KEY })
}

/**
 * Refreshes related queries after creating an experiment.
 *
 * Input example:
 * - `{ id: 21, name: 'Dose response', project: 5 }`
 */
const onExperimentCreated = async (_experiment: Experiment) => {
  await Promise.all([
    queryClient.invalidateQueries({ queryKey: PROJECTS_QUERY_KEY }),
    queryClient.invalidateQueries({ queryKey: EXPERIMENTS_QUERY_KEY }),
  ])
}
</script>

<template>
  <section class="space-y-1">
    <p class="px-2 text-xs font-semibold tracking-wide text-slate-500 uppercase">
      {{ t('navigation.projects.title') }}
    </p>
    <p v-if="projectsQuery.isPending.value" class="px-2 text-xs text-slate-500">
      {{ t('navigation.projects.loading') }}
    </p>
    <p v-else-if="projectsErrorMessage" class="px-2 text-xs text-red-600">
      {{ projectsErrorMessage }}
    </p>
    <p v-else-if="projectNodes.length === 0" class="px-2 text-xs text-slate-500">
      {{ t('navigation.projects.empty') }}
    </p>

    <UTree v-else :items="items" :ui="{ link: 'cursor-pointer hover:text-primary transition-colors' }">
      <template #projects-root-label="{ item }">
        <UContextMenu :items="projectsRootContextMenuItems">
          <span class="inline-block w-full">{{ getSlotItemLabel(item) }}</span>
        </UContextMenu>
      </template>

      <template #project-node-label="{ item }">
        <UContextMenu :items="getProjectContextMenuItemsFromSlotItem(item)">
          <span class="inline-block w-full">{{ getSlotItemLabel(item) }}</span>
        </UContextMenu>
      </template>
    </UTree>

    <ProjectCreateModal v-model:open="isCreateProjectModalOpen" @created="onProjectCreated" />
    <ExperimentCreateModal
      v-model:open="isCreateExperimentModalOpen"
      :project-id="selectedProjectIdForNewExperiment"
      :project-name="selectedProjectNameForNewExperiment"
      @created="onExperimentCreated"
    />
  </section>
</template>
