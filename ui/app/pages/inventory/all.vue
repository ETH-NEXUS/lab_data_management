<script setup lang="ts">
import { computed } from 'vue'
import InventoryStockWorkspace from '~/components/inventory/InventoryStockWorkspace.vue'

type StockPreset = 'all' | 'favorite' | 'low_stock' | 'expired'

const { t } = useI18n()
const route = useRoute()

const goBackToInventory = (): void => {
  navigateTo('/inventory')
}

/**
 * Normalizes `preset` query values to one allowed stock workspace preset.
 *
 * Accepted query examples:
 * - `?preset=all`
 * - `?preset=favorite`
 *
 * Returned examples:
 * - `'all'`
 * - `'favorite'`
 * - fallback `'all'` for unknown/missing values
 */
const selectedPreset = computed<StockPreset>(() => {
  const presetValue = route.query.preset
  const normalizedPreset = Array.isArray(presetValue) ? presetValue[0] : presetValue

  if (normalizedPreset === 'favorite') return 'favorite'
  if (normalizedPreset === 'low_stock') return 'low_stock'
  if (normalizedPreset === 'expired') return 'expired'
  return 'all'
})
</script>

<template>
  <section class="inventory-shell px-2 pb-12 sm:px-3 lg:px-4">
    <div class="w-full space-y-5">
      <UButton
        variant="ghost"
        color="neutral"
        icon="i-heroicons-arrow-left"
        :label="t('inventory.page.title')"
        @click="goBackToInventory"
      />

      <InventoryStockWorkspace :preset="selectedPreset" />
    </div>
  </section>
</template>
