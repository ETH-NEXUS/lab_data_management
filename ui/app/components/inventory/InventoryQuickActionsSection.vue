<script setup lang="ts">
import { computed } from 'vue'

type InventoryActionItem = {
  id: string
  title: string
  icon?: string
}

const emit = defineEmits<{
  (e: 'select-action', actionId: string): void
}>()

const { t } = useI18n()

/**
 * Primary "all items" action rendered as direct link button.
 *
 * Returned data example:
 * - `{ id: 'all_items', title: 'All items' }`
 */
const allItemsAction = computed<InventoryActionItem>(() => ({
  id: 'all_items',
  title: t('inventory.page.actions.all_items.title'),
}))

/**
 * Top-row operational actions (compact card style).
 *
 * Returned data example:
 * - `[{ id: 'search', title: 'Search', icon: 'i-heroicons-magnifying-glass' }]`
 */
const topRowActions = computed<InventoryActionItem[]>(() => [
  {
    id: 'search',
    title: t('inventory.page.actions.search.title'),
    icon: 'i-heroicons-magnifying-glass',
  },
  {
    id: 'advanced_search',
    title: t('inventory.page.actions.advanced_search.title'),
    icon: 'i-heroicons-magnifying-glass',
  },
  {
    id: 'order',
    title: t('inventory.page.actions.order.title'),
  },
])

/**
 * Secondary icon actions shown separately to the right.
 *
 * Returned data example:
 * - `[{ id: 'add_new_item', title: 'Add new item', icon: 'i-heroicons-plus' }]`
 */
const sideActions = computed<InventoryActionItem[]>(() => [
  {
    id: 'add_new_item',
    title: t('inventory.page.actions.add_new_item.title'),
    icon: 'i-heroicons-plus',
  },
  {
    id: 'favorite_items',
    title: t('inventory.page.actions.favorite_items.title'),
    icon: 'i-heroicons-star',
  },
])

/**
 * Forwards selected quick action id to parent layout.
 *
 * Input example:
 * - `actionId = 'all_items'`
 */
const onSelectAction = (actionId: string): void => {
  emit('select-action', actionId)
}
</script>

<template>
  <section class="space-y-4">
    <h2 class="inventory-compact-title">{{ t('inventory.page.title') }}</h2>

    <div class="grid gap-4 lg:grid-cols-[minmax(0,1fr)_auto] lg:items-start">
      <div class="grid gap-3 sm:grid-cols-2 xl:grid-cols-4">
        <NuxtLink to="/inventory/all" class="inventory-toolbar-card inventory-toolbar-card--primary">
          <span class="truncate">{{ allItemsAction.title }}</span>
        </NuxtLink>

        <button
          v-for="action in topRowActions"
          :key="action.id"
          type="button"
          class="inventory-toolbar-card"
          @click="onSelectAction(action.id)"
        >
          <span class="truncate">{{ action.title }}</span>
          <UIcon v-if="action.icon" :name="action.icon" class="size-5 shrink-0 text-slate-600" />
        </button>
      </div>

      <div class="flex gap-3 lg:flex-col">
        <button
          v-for="action in sideActions"
          :key="action.id"
          type="button"
          class="inventory-side-card"
          @click="onSelectAction(action.id)"
        >
          <span class="inventory-side-card-icon">
            <UIcon :name="action.icon || 'i-heroicons-circle-stack'" class="size-5 text-slate-700" />
          </span>
          <span class="text-sm text-slate-700">{{ action.title }}</span>
        </button>
      </div>
    </div>
  </section>
</template>
