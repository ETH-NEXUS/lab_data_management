<script setup lang="ts">
import { computed } from 'vue'
import InventoryActionCard from '~/components/inventory/InventoryActionCard.vue'

type InventoryActionItem = {
  id: string
  title: string
  description: string
  icon: string
}

const emit = defineEmits<{
  (e: 'select-action', actionId: string): void
}>()

const { t } = useI18n()

/**
 * Creates quick-action cards shown in the first inventory action row.
 *
 * Returned data example:
 * - `[{ id: 'search', title: 'Search', description: 'Find items...', icon: 'i-heroicons-magnifying-glass' }]`
 */
const quickActions = computed<InventoryActionItem[]>(() => [
  {
    id: 'all_items',
    title: t('inventory.page.actions.all_items.title'),
    description: t('inventory.page.actions.all_items.description'),
    icon: 'i-heroicons-squares-2x2',
  },
  {
    id: 'search',
    title: t('inventory.page.actions.search.title'),
    description: t('inventory.page.actions.search.description'),
    icon: 'i-heroicons-magnifying-glass',
  },
  {
    id: 'advanced_search',
    title: t('inventory.page.actions.advanced_search.title'),
    description: t('inventory.page.actions.advanced_search.description'),
    icon: 'i-heroicons-adjustments-horizontal',
  },
  {
    id: 'order',
    title: t('inventory.page.actions.order.title'),
    description: t('inventory.page.actions.order.description'),
    icon: 'i-heroicons-shopping-cart',
  },
  {
    id: 'add_new_item',
    title: t('inventory.page.actions.add_new_item.title'),
    description: t('inventory.page.actions.add_new_item.description'),
    icon: 'i-heroicons-plus-circle',
  },
  {
    id: 'favorite_items',
    title: t('inventory.page.actions.favorite_items.title'),
    description: t('inventory.page.actions.favorite_items.description'),
    icon: 'i-heroicons-star',
  },
])

/**
 * Forwards a selected quick-action id to the parent inventory layout.
 *
 * Input example:
 * - `actionId = 'order'`
 */
const onSelectAction = (actionId: string): void => {
  emit('select-action', actionId)
}
</script>

<template>
  <section class="space-y-3">
    <p class="inventory-section-title">{{ t('inventory.page.section_labels.quick_actions') }}</p>

    <div class="grid gap-4 sm:grid-cols-2 xl:grid-cols-3">
      <InventoryActionCard
        v-for="action in quickActions"
        :key="action.id"
        :item="action"
        size="quick"
        @select="onSelectAction"
      />
    </div>
  </section>
</template>
