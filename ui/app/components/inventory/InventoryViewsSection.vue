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
 * Creates inventory view cards shown in the main inventory grid.
 *
 * Returned data example:
 * - `[{ id: 'devices', title: 'Devices', description: 'Manage equipment...', icon: 'i-heroicons-computer-desktop' }]`
 */
const inventoryViews = computed<InventoryActionItem[]>(() => [
  {
    id: 'expired_items',
    title: t('inventory.page.actions.expired_items.title'),
    description: t('inventory.page.actions.expired_items.description'),
    icon: 'i-heroicons-clock',
  },
  {
    id: 'devices',
    title: t('inventory.page.actions.devices.title'),
    description: t('inventory.page.actions.devices.description'),
    icon: 'i-heroicons-computer-desktop',
  },
  {
    id: 'recent_activities',
    title: t('inventory.page.actions.recent_activities.title'),
    description: t('inventory.page.actions.recent_activities.description'),
    icon: 'i-heroicons-bolt',
  },
  {
    id: 'recently_linked_items',
    title: t('inventory.page.actions.recently_linked_items.title'),
    description: t('inventory.page.actions.recently_linked_items.description'),
    icon: 'i-heroicons-link',
  },
  {
    id: 'reserved_items',
    title: t('inventory.page.actions.reserved_items.title'),
    description: t('inventory.page.actions.reserved_items.description'),
    icon: 'i-heroicons-bookmark',
  },
  {
    id: 'archived_items',
    title: t('inventory.page.actions.archived_items.title'),
    description: t('inventory.page.actions.archived_items.description'),
    icon: 'i-heroicons-archive-box',
  },
])

/**
 * Forwards one selected inventory view id to the parent inventory layout.
 *
 * Input example:
 * - `actionId = 'archived_items'`
 */
const onSelectAction = (actionId: string): void => {
  emit('select-action', actionId)
}
</script>

<template>
  <section class="space-y-3">
    <p class="inventory-section-title">{{ t('inventory.page.section_labels.inventory_views') }}</p>

    <div class="grid gap-4 sm:grid-cols-2 xl:grid-cols-3">
      <InventoryActionCard
        v-for="action in inventoryViews"
        :key="action.id"
        :item="action"
        size="view"
        @select="onSelectAction"
      />
    </div>
  </section>
</template>
