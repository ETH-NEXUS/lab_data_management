<script setup lang="ts">
import { computed } from 'vue'
import type { NavigationTreeNode } from '~/types/navigation'

const props = withDefaults(defineProps<{ filter?: string }>(), {
  filter: '',
})

const rawItems = computed<NavigationTreeNode[]>(() => [
  {
    label: 'Inventory',
    icon: 'i-heroicons-cube',
    onSelect: () => navigateTo('/inventory'),
  },
])

/**
 * Filters inventory tree items by a case-insensitive label match.
 *
 * Input example:
 * - `query = 'inv'`
 *
 * Returned example:
 * - `[{ label: 'Inventory', icon: 'i-heroicons-cube', onSelect: () => void }]`
 */
const items = computed<NavigationTreeNode[]>(() => {
  const query = props.filter.trim().toLowerCase()
  if (query === '') {
    return rawItems.value
  }

  const filteredItems: NavigationTreeNode[] = []
  for (const item of rawItems.value) {
    if (item.label.toLowerCase().includes(query)) {
      filteredItems.push(item)
    }
  }

  return filteredItems
})
</script>

<template>
  <section v-if="items.length > 0" class="space-y-1">
    <p class="px-2 text-xs font-semibold tracking-wide text-slate-500 uppercase">Inventory</p>
    <UTree :items="items" :ui="{ link: 'cursor-pointer hover:text-blue-700 transition-colors' }" />
  </section>
</template>
