<script setup lang="ts">
import { computed } from 'vue'
import InventoryStockWorkspace from '~/components/inventory/InventoryStockWorkspace.vue'
import type { InventoryStockPreset } from '~/types/inventory'

const { t } = useI18n()
const route = useRoute()
const router = useRouter()

const goBackToInventory = (): void => {
  navigateTo('/inventory')
}

/**
 * Normalizes `preset` query values to one allowed stock workspace preset.
 *
 * Accepted query examples:
 * - `?preset=all`
 * - `?preset=favorite`
 * - `?preset=archived`
 *
 * Returned examples:
 * - `'all'`
 * - `'favorite'`
 * - `'archived'`
 * - fallback `'all'` for unknown/missing values
 */
const selectedPreset = computed<InventoryStockPreset | undefined>(() => {
  const presetValue = route.query.preset
  const normalizedPreset = Array.isArray(presetValue) ? presetValue[0] : presetValue

  if (normalizedPreset === 'favorite') return 'favorite'
  if (normalizedPreset === 'low_stock') return 'low_stock'
  if (normalizedPreset === 'expired') return 'expired'
  if (normalizedPreset === 'archived') return 'archived'
  if (normalizedPreset === 'all') return 'all'
  return undefined
})

const selectedStockId = computed<number | null>(() => {
  const stockValue = route.query.stock
  const normalizedStockValue = Array.isArray(stockValue) ? stockValue[0] : stockValue

  if (!normalizedStockValue) {
    return null
  }

  const parsedStockId = Number.parseInt(normalizedStockValue, 10)

  return Number.isNaN(parsedStockId) ? null : parsedStockId
})

/**
 * Keeps the `stock` query in sync with the selected stock drawer state.
 *
 * Input examples:
 * - `14`
 * - `null`
 */
const updateSelectedStockId = async (stockId: number | null): Promise<void> => {
  const nextQuery = { ...route.query }

  if (stockId === null) {
    delete nextQuery.stock
  } else {
    nextQuery.stock = String(stockId)
  }

  await router.replace({
    query: nextQuery,
  })
}
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

      <InventoryStockWorkspace
        :preset="selectedPreset"
        :initial-stock-id="selectedStockId"
        @update:selected-stock-id="updateSelectedStockId"
      />
    </div>
  </section>
</template>
