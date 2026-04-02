<script setup lang="ts">
import InventoryLandingContent from '~/components/inventory/InventoryLandingContent.vue'

const { t } = useI18n()

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
  const titleKey = `inventory.page.actions.${actionId}.title`
  const actionTitle = t(titleKey)
  const placeholderSuffix = t('inventory.page.placeholder_suffix')

  window.alert(`${actionTitle}: ${placeholderSuffix}`)
}
</script>

<template>
  <InventoryLandingContent @select-action="onSelectAction" />
</template>
