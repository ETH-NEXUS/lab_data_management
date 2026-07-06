<script setup lang="ts">
import { ref } from 'vue'
import InventoryAddItemModal from '~/components/inventory/InventoryAddItemModal.vue'
import InventoryLandingContent from '~/components/inventory/InventoryLandingContent.vue'

const { t } = useI18n()
const isAddItemModalOpen = ref(false)

const onUpdateAddItemModalOpen = (isOpen: boolean): void => {
  isAddItemModalOpen.value = isOpen
}

/**
 * Shows a temporary placeholder alert for selected inventory actions.
 *
 * Accepted data example:
 * - `actionId = 'devices'`
 *
 * Returned behavior example:
 * - `window.alert('Devices: to be implemented')`
 */
const onSelectAction = (actionId: string): void => {
  if (actionId === 'all_items') {
    navigateTo('/inventory/all?preset=all')
    return
  }

  if (actionId === 'add_new_item') {
    isAddItemModalOpen.value = true
    return
  }

  if (actionId === 'recently_linked_items') {
    navigateTo('/inventory/usages')
    return
  }

  if (actionId === 'favorite_items') {
    navigateTo('/inventory/all?preset=favorite')
    return
  }

  if (actionId === 'expired_items') {
    navigateTo('/inventory/all?preset=expired')
    return
  }

  if (actionId === 'low_stock_items') {
    navigateTo('/inventory/all?preset=low_stock')
    return
  }

  if (actionId === 'archived_items') {
    navigateTo('/inventory/all?preset=archived')
    return
  }

  const titleKey = `inventory.page.actions.${actionId}.title`
  const actionTitle = t(titleKey)
  const placeholderSuffix = t('inventory.page.placeholder_suffix')

  window.alert(`${actionTitle}: ${placeholderSuffix}`)
}
</script>

<template>
  <div>
    <InventoryLandingContent @select-action="onSelectAction" />
    <InventoryAddItemModal :open="isAddItemModalOpen" @update:open="onUpdateAddItemModalOpen" />
  </div>
</template>
