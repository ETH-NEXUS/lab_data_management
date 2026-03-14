<script setup lang="ts">
import { computed } from 'vue'

type TreeNode = {
  label: string
  icon?: string
  defaultExpanded?: boolean
  children?: TreeNode[]
  onSelect?: () => void
}

const props = withDefaults(defineProps<{ filter?: string }>(), {
  filter: '',
})

const rawItems = computed<TreeNode[]>(() => [
  {
    label: 'Management',
    icon: 'i-heroicons-cog-6-tooth',
    defaultExpanded: true,
    children: [
      {
        label: 'Users',
        icon: 'i-heroicons-users',
        onSelect: () => navigateTo('/management/users'),
      },
      {
        label: 'System Logs',
        icon: 'i-heroicons-clipboard-document-list',
        onSelect: () => navigateTo('/management/logs'),
      },
      {
        label: 'Task Runner',
        icon: 'i-heroicons-play',
        onSelect: () => navigateTo('/management/tasks'),
      },
    ],
  },
])

const filterTree = (items: TreeNode[], query: string): TreeNode[] => {
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
    .filter((item): item is TreeNode => item !== null)
}

const items = computed<TreeNode[]>(() => filterTree(rawItems.value, props.filter))
</script>

<template>
  <section class="space-y-1">
    <p class="px-2 text-xs font-semibold uppercase tracking-wide text-slate-500">
      Management
    </p>
    <UTree
      :items="items"
      :ui="{ link: 'cursor-pointer hover:text-primary transition-colors' }"
    />
  </section>
</template>
