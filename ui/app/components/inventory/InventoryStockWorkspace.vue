<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import InventoryStockDetailDrawer from '~/components/inventory/InventoryStockDetailDrawer.vue'
import InventoryStockTable from '~/components/inventory/InventoryStockTable.vue'
import { useInventoryMaterialQuery } from '~/composables/inventory/useInventoryMaterialQuery'
import { useInventoryOrdersQuery } from '~/composables/inventory/useInventoryOrderQuery'
import { useInventoryStocksQuery } from '~/composables/inventory/useInventoryStockQuery'
import { useInventoryUsagesQuery } from '~/composables/inventory/useInventoryUsageQuery'
import type { InventoryOrderListItem, InventoryStockListItem, InventoryUsageListItem } from '~/types/inventory'

const { t } = useI18n()

const stocksQuery = useInventoryStocksQuery()
const usagesQuery = useInventoryUsagesQuery()
const ordersQuery = useInventoryOrdersQuery()

const selectedStockId = ref<number | null>(null)

const stocks = computed<InventoryStockListItem[]>(() => stocksQuery.data.value ?? [])

const selectedStock = computed<InventoryStockListItem | null>(() => {
  if (!selectedStockId.value) {
    return null
  }

  for (const stock of stocks.value) {
    if (stock.id === selectedStockId.value) {
      return stock
    }
  }

  return null
})

const selectedMaterialId = computed<number>(() => {
  if (!selectedStock.value) {
    return 0
  }
  return selectedStock.value.material.id
})

const selectedMaterialQuery = useInventoryMaterialQuery(selectedMaterialId)

const selectedStockUsageEntries = computed<InventoryUsageListItem[]>(() => {
  const stock = selectedStock.value
  if (!stock) {
    return []
  }

  const usages = usagesQuery.data.value ?? []
  const matches: InventoryUsageListItem[] = []

  for (const usage of usages) {
    if (usage.inventory_stock.id === stock.id) {
      matches.push(usage)
    }
  }

  return matches
})

const selectedStockOrderEntries = computed<InventoryOrderListItem[]>(() => {
  const stock = selectedStock.value
  if (!stock) {
    return []
  }

  const orders = ordersQuery.data.value ?? []
  const matches: InventoryOrderListItem[] = []

  for (const order of orders) {
    if (order.material.id === stock.material.id) {
      matches.push(order)
    }
  }

  return matches
})

const stocksErrorMessage = computed<string | null>(() => {
  const err = stocksQuery.error.value
  if (!err) {
    return null
  }
  return err instanceof Error ? err.message : t('inventory.stock_workspace.error')
})

/**
 * Opens the stock detail drawer for one selected row.
 *
 * Input example:
 * - `{ id: 14, material: { id: 3, product_name: 'PBS buffer' }, ... }`
 */
const openStockDetails = (stock: InventoryStockListItem): void => {
  selectedStockId.value = stock.id
}

const closeStockDetails = (): void => {
  selectedStockId.value = null
}

// Close drawer when selected stock no longer exists in refreshed data.
watch(
  () => selectedStockId.value,
  (stockId) => {
    if (!stockId) {
      return
    }

    const hasMatch = stocks.value.some((stock) => stock.id === stockId)
    if (!hasMatch) {
      selectedStockId.value = null
    }
  },
)

watch(
  () => selectedStock.value,
  (stock) => {
    if (!stock && selectedStockId.value !== null) {
      selectedStockId.value = null
    }
  },
)
</script>

<template>
  <section class="space-y-3">
    <p class="inventory-section-title">{{ t('inventory.stock_workspace.title') }}</p>

    <UCard
      :ui="{
        root: 'core-card divide-y divide-slate-200/70',
      }"
    >
      <template #header>
        <div class="space-y-1">
          <p class="text-sm font-semibold text-slate-800">{{ t('inventory.stock_workspace.title') }}</p>
          <p class="text-sm text-slate-600">{{ t('inventory.stock_workspace.description') }}</p>
        </div>
      </template>

      <p v-if="stocksQuery.isPending.value" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.loading') }}
      </p>
      <p v-else-if="stocksErrorMessage" class="text-sm text-red-600">
        {{ stocksErrorMessage }}
      </p>
      <p v-else-if="stocks.length === 0" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.empty') }}
      </p>
      <InventoryStockTable v-else :stocks="stocks" @select-stock="openStockDetails" />
    </UCard>

    <InventoryStockDetailDrawer
      :open="selectedStock !== null"
      :stock="selectedStock"
      :material-detail="selectedMaterialQuery.data.value ?? null"
      :usage-entries="selectedStockUsageEntries"
      :order-entries="selectedStockOrderEntries"
      :is-material-loading="selectedMaterialQuery.isPending.value"
      :is-usages-loading="usagesQuery.isPending.value"
      :is-orders-loading="ordersQuery.isPending.value"
      @close="closeStockDetails"
    />
  </section>
</template>
