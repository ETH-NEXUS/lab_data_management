<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import InventoryStockDetailsPanel from '~/components/inventory/stock-details/InventoryStockDetailsPanel.vue'
import InventoryStockTable from '~/components/inventory/InventoryStockTable.vue'
import { useInventoryStocksQuery } from '~/composables/inventory/useInventoryStockQuery'
import type { InventoryStockListItem } from '~/types/inventory'

type StockPreset = 'all' | 'favorite' | 'low_stock' | 'expired'

type Props = {
  preset?: StockPreset
}

const props = withDefaults(defineProps<Props>(), {
  preset: 'all',
})

const { t } = useI18n()

const stocksQuery = useInventoryStocksQuery()

const selectedStockId = ref<number | null>(null)

const stocks = computed<InventoryStockListItem[]>(() => stocksQuery.data.value ?? [])

/**
 * Applies one preset filter on top of the already-loaded full stocks list.
 *
 * Input examples:
 * - `preset = 'favorite'`
 * - `preset = 'low_stock'`
 *
 * Returned data example:
 * - `[{ id: 4, is_favorite: true, ... }, { id: 12, is_favorite: true, ... }]`
 */
const filteredStocks = computed<InventoryStockListItem[]>(() => {
  if (props.preset === 'favorite') {
    return stocks.value.filter((stock) => stock.is_favorite)
  }

  if (props.preset === 'low_stock') {
    return stocks.value.filter((stock) => stock.is_low_stock || stock.inventory_status === 'low')
  }

  if (props.preset === 'expired') {
    return stocks.value.filter((stock) => stock.is_expired)
  }

  return stocks.value
})

const selectedStock = computed<InventoryStockListItem | null>(() => {
  if (!selectedStockId.value) {
    return null
  }

  for (const stock of filteredStocks.value) {
    if (stock.id === selectedStockId.value) {
      return stock
    }
  }

  return null
})

const stocksErrorMessage = computed<string | null>(() => {
  const err = stocksQuery.error.value
  if (!err) {
    return null
  }
  return err instanceof Error ? err.message : t('inventory.stock_workspace.error')
})

const hasSelectedStock = computed<boolean>(() => {
  return selectedStock.value !== null
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

    const hasMatch = filteredStocks.value.some((stock) => stock.id === stockId)
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

    <div
      :class="[
        'grid items-start gap-4',
        hasSelectedStock ? 'xl:grid-cols-[minmax(0,1.45fr)_minmax(0,1fr)]' : 'grid-cols-1',
      ]"
    >
      <UCard
        class="min-w-0"
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
        <p v-else-if="filteredStocks.length === 0" class="text-sm text-slate-600">
          {{ t('inventory.stock_workspace.empty') }}
        </p>
        <InventoryStockTable v-else :stocks="filteredStocks" @select-stock="openStockDetails" />
      </UCard>

      <InventoryStockDetailsPanel
        v-if="hasSelectedStock"
        class="min-w-0"
        :open="hasSelectedStock"
        :stock="selectedStock"
        @close="closeStockDetails"
      />
    </div>
  </section>
</template>
