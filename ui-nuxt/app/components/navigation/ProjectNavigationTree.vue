<script setup lang="ts">
import { computed } from 'vue'
import { useProjectsQuery } from '~/composables/useProjectsQuery'
import type { NavigationTreeNode } from '~/types/navigation'
import { PROJECTS_ERROR_MESSAGE } from '~/types/projects'

const props = withDefaults(defineProps<{ filter?: string }>(), {
  filter: '',
})

const projectsQuery = useProjectsQuery()

const projectNodes = computed<NavigationTreeNode[]>(() =>
  (projectsQuery.data.value ?? []).map((project) => ({
    label: project.name || `Project ${project.id}`,
    icon: 'i-heroicons-beaker',
  })),
)

const rawItems = computed<NavigationTreeNode[]>(() => [
  {
    label: 'Projects',
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
  return err instanceof Error ? err.message : PROJECTS_ERROR_MESSAGE
})
</script>

<template>
  <section class="space-y-1">
    <p class="px-2 text-xs font-semibold tracking-wide text-slate-500 uppercase">Projects</p>

    <p v-if="projectsQuery.isPending.value" class="px-2 text-xs text-slate-500">Loading projects...</p>
    <p v-else-if="projectsErrorMessage" class="px-2 text-xs text-red-600">
      {{ projectsErrorMessage }}
    </p>
    <p v-else-if="projectNodes.length === 0" class="px-2 text-xs text-slate-500">No projects found.</p>
    <UTree v-else :items="items" :ui="{ link: 'cursor-pointer hover:text-primary transition-colors' }" />
  </section>
</template>
