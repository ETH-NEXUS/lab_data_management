<script setup lang="ts">
import InventoryQuickActionsSection from '~/components/inventory/InventoryQuickActionsSection.vue'
import InventoryViewsSection from '~/components/inventory/InventoryViewsSection.vue'

const emit = defineEmits<{
  (e: 'select-action', actionId: string): void
}>()

/**
 * Passes selected action ids from section components to the page level.
 *
 * Input example:
 * - `actionId = 'favorite_items'`
 *
 * Returned data example:
 * - emits `select-action` with `'favorite_items'`
 */
const onSelectAction = (actionId: string): void => {
  emit('select-action', actionId)
}
</script>

<template>
  <section class="inventory-shell px-4 sm:px-6 lg:px-8">
    <div class="mx-auto w-full max-w-6xl space-y-7">
      <InventoryQuickActionsSection @select-action="onSelectAction" />
      <InventoryViewsSection @select-action="onSelectAction" />
    </div>
  </section>
</template>
