<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import InventoryStockDetailsPanel from '~/components/inventory/stock-details/InventoryStockDetailsPanel.vue'
import InventoryStockTable from '~/components/inventory/InventoryStockTable.vue'
import { getStocksForPreset } from '~/components/inventory/inventory-stock-table.values'
import { useInventoryStocksQuery } from '~/composables/inventory/useInventoryStockQuery'
import { useInventoryStockTablePreferenceStore } from '~/stores/inventory/InventoryStockTablePreferenceStore'
import type { InventoryStockListItem, InventoryStockPreset } from '~/types/inventory'

type Props = {
  preset?: InventoryStockPreset
  initialStockId?: number | null
}

const props = defineProps<Props>()
const emit = defineEmits<{
  (e: 'update:selected-stock-id', stockId: number | null): void
}>()

const { t } = useI18n()

const stocksQuery = useInventoryStocksQuery()
const stockTablePreferenceStore = useInventoryStockTablePreferenceStore()

const selectedStockId = ref<number | null>(props.initialStockId ?? null)

const stocks = computed<InventoryStockListItem[]>(() => stocksQuery.data.value ?? [])
const effectivePreset = computed<InventoryStockPreset>(() => {
  return props.preset ?? stockTablePreferenceStore.presetState
})

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
  return getStocksForPreset(stocks.value, effectivePreset.value)
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

const workspaceTitle = computed<string>(() => {
  if (effectivePreset.value === 'favorite') {
    return t('inventory.page.actions.favorite_items.title')
  }

  if (effectivePreset.value === 'low_stock') {
    return t('inventory.page.actions.low_stock_items.title')
  }

  if (effectivePreset.value === 'expired') {
    return t('inventory.page.actions.expired_items.title')
  }

  return t('inventory.stock_workspace.title')
})

const workspaceDescription = computed<string>(() => {
  if (effectivePreset.value === 'favorite') {
    return t('inventory.page.actions.favorite_items.description')
  }

  if (effectivePreset.value === 'low_stock') {
    return t('inventory.page.actions.low_stock_items.description')
  }

  if (effectivePreset.value === 'expired') {
    return t('inventory.page.actions.expired_items.description')
  }

  return t('inventory.stock_workspace.description')
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
  emit('update:selected-stock-id', stock.id)
}

const closeStockDetails = (): void => {
  selectedStockId.value = null
  emit('update:selected-stock-id', null)
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
      emit('update:selected-stock-id', null)
    }
  },
)

watch(
  () => selectedStock.value,
  (stock) => {
    if (!stock && selectedStockId.value !== null) {
      selectedStockId.value = null
      emit('update:selected-stock-id', null)
    }
  },
)

watch(
  () => props.preset,
  (nextPreset) => {
    if (!nextPreset) {
      return
    }

    if (nextPreset !== stockTablePreferenceStore.presetState) {
      void stockTablePreferenceStore.updatePresetState(nextPreset)
    }
  },
  { immediate: true },
)

watch(
  () => props.initialStockId,
  (nextStockId) => {
    selectedStockId.value = nextStockId ?? null
  },
  { immediate: true },
)
</script>

<template>
  <section class="space-y-3">
    <p class="inventory-section-title">{{ workspaceTitle }}</p>

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
            <p class="text-sm font-semibold text-slate-800">{{ workspaceTitle }}</p>
            <p class="text-sm text-slate-600">{{ workspaceDescription }}</p>
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
        <InventoryStockTable
          v-else
          :stocks="filteredStocks"
          :sorting-state="stockTablePreferenceStore.sortingState"
          :global-filter-state="stockTablePreferenceStore.globalFilterState"
          :column-filters-state="stockTablePreferenceStore.columnFiltersState"
          :column-order-state="stockTablePreferenceStore.columnOrderState"
          :column-visibility-state="stockTablePreferenceStore.columnVisibilityState"
          @sorting-change="stockTablePreferenceStore.updateSortingState"
          @global-filter-change="stockTablePreferenceStore.updateGlobalFilterState"
          @column-filters-change="stockTablePreferenceStore.updateColumnFiltersState"
          @column-order-change="stockTablePreferenceStore.updateColumnOrderState"
          @column-visibility-change="stockTablePreferenceStore.updateColumnVisibilityState"
          @select-stock="openStockDetails"
        />
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
