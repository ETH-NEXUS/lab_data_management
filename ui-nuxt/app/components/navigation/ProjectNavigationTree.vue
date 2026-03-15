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
/**
 * Creates plate nodes for one experiment.
 *
 * Accepted data example:
 * - `experiment = { id: 21, name: 'Dose response', plates: [{ id: 11, barcode: 'A001', dimension: { name: '384' } }] }`
 *
 * Returned data example:
 * - `[{ label: 'A001 (384)', icon: 'i-heroicons-squares-2x2', onSelect: () => void }]`
 */
const mapPlateNodes = (experiment: Experiment): NavigationTreeNode[] => {
  // Step 1: build a plain array (avoid compact spread syntax for readability).
  const plates = experiment.plates ? [...experiment.plates] : []

  // Step 2: sort by barcode so users see stable ordering.
  plates.sort((leftPlate, rightPlate) => leftPlate.barcode.localeCompare(rightPlate.barcode))

  // Step 3: map each plate to a tree node.
  const nodes: NavigationTreeNode[] = []
  for (const plate of plates) {
    const plateLabel = plate.barcode || t('navigation.projects.plate_fallback', { id: plate.id })
    const dimensionLabel = formatPlateDimensionLabel(plate.dimension)
    const routePlateId = plate.barcode || String(plate.id)

    nodes.push({
      label: `${plateLabel} (${dimensionLabel})`,
      icon: 'i-heroicons-squares-2x2',
      onSelect: () => navigateTo(`/plates/${encodeURIComponent(routePlateId)}`),
    })
  }

  return nodes
}

/**
 * Creates control-plate nodes for one project.
 *
 * Accepted data example:
 * - `project = { id: 5, name: 'Screening', plates: [{ id: 2, barcode: 'CTRL-001' }] }`
 *
 * Returned data example:
 * - `[{ label: 'Control plate: CTRL-001', icon: 'i-heroicons-squares-2x2', onSelect: () => void }]`
 */
const mapControlPlateNodes = (project: Project): NavigationTreeNode[] => {
  // make an explicit array from optional project.plates.
  const controlPlates = project.plates ? [...project.plates] : []

  // sort by barcode for predictable order in the tree.
  controlPlates.sort((leftPlate, rightPlate) => leftPlate.barcode.localeCompare(rightPlate.barcode))

  // : create tree nodes.
  const nodes: NavigationTreeNode[] = []
  for (const plate of controlPlates) {
    const plateLabel = plate.barcode || t('navigation.projects.plate_fallback', { id: plate.id })
    const routePlateId = plate.barcode || String(plate.id)

    nodes.push({
      label: t('navigation.projects.control_plate', { barcode: plateLabel }),
      icon: 'i-heroicons-squares-2x2',
      onSelect: () => navigateTo(`/plates/${encodeURIComponent(routePlateId)}`),
    })
  }

  return nodes
}

/**
 * Creates experiment nodes (with plate children) for one project.
 *
 * Accepted data example:
 * - `project = { id: 5, experiments: [{ id: 21, name: 'Dose response', plates: [] }] }`
 *
 * Returned data example:
 * - `[{ label: 'Dose response', icon: 'i-heroicons-beaker', children: [] }]`
 */
const mapExperimentNodes = (project: Project): NavigationTreeNode[] => {
  const experiments = project.experiments ?? []
  const nodes: NavigationTreeNode[] = []

  for (const experiment of experiments) {
    const experimentLabel = experiment.name || t('navigation.projects.experiment_fallback', { id: experiment.id })

    nodes.push({
      label: experimentLabel,
      icon: 'i-heroicons-beaker',
      onSelect: () => navigateTo(`/experiments/${experiment.id}`),
      children: mapPlateNodes(experiment),
    })
  }

  return nodes
}

/**
 * Builds all project children in the same order as legacy UI:
 * 1. control plates
 * 2. experiments
 *
 * Returned data example:
 * - `[{ label: 'Control plate: CTRL-001' }, { label: 'Dose response' }]`
 */
const mapProjectChildNodes = (project: Project): NavigationTreeNode[] => {
  const nodes: NavigationTreeNode[] = []

  const controlPlateNodes = mapControlPlateNodes(project)
  for (const node of controlPlateNodes) {
    nodes.push(node)
  }

  const experimentNodes = mapExperimentNodes(project)
  for (const node of experimentNodes) {
    nodes.push(node)
  }

  return nodes
}

/**
 * Creates first-level project nodes for UTree.
 *
 * Returned data example:
 * - `[{ label: 'Screening 2026', projectId: 5, children: [...] }]`
 */
const projectNodes = computed<NavigationTreeNode[]>(() => {
  const projects = projectsQuery.data.value ?? []
  const nodes: NavigationTreeNode[] = []

  for (const project of projects) {
    const projectLabel = project.name || t('navigation.projects.project_fallback', { id: project.id })

    nodes.push({
      label: projectLabel,
      slot: 'project-node',
      projectId: project.id,
      projectName: projectLabel,
      icon: 'i-heroicons-document-text',
      onSelect: () => navigateTo(`/projects/${project.id}`),
      children: mapProjectChildNodes(project),
    })
  }

  return nodes
})

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

  const filteredItems: NavigationTreeNode[] = []

  for (const item of items) {
    const filteredChildren = item.children ? filterTree(item.children, q) : []
    const selfMatches = item.label.toLowerCase().includes(q)
    const hasMatchingChildren = filteredChildren.length > 0

    if (!selfMatches && !hasMatchingChildren) {
      continue
    }

    filteredItems.push({
      ...item,
      children: filteredChildren,
      defaultExpanded: true,
    })
  }

  return filteredItems
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
  [
    {
      label: t('navigation.projects.actions.select_control_layout'),
      icon: 'i-heroicons-squares-2x2',
      onSelect: () => navigateTo(`/layout/${projectId}`),
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
