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
    label: 'Templates',
    icon: 'i-heroicons-document-duplicate',
    defaultExpanded: true,
    children: [
      {
        label: 'Plate Layout Templates',
        icon: 'i-heroicons-table-cells',
        onSelect: () => navigateTo('/templates/plate-layouts'),
      },
      {
        label: 'Report Templates',
        icon: 'i-heroicons-document-text',
        onSelect: () => navigateTo('/templates/reports'),
      },
      {
        label: 'Legacy Template Import',
        icon: 'i-heroicons-arrow-up-tray',
        onSelect: () => navigateTo('/templates/import'),
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
      Templates
    </p>
    <UTree
      :items="items"
      :ui="{ link: 'cursor-pointer hover:text-primary transition-colors' }"
    />
  </section>
</template>
