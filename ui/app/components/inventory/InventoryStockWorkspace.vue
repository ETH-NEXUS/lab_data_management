<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import type { PaginationState, SortingState } from '@tanstack/vue-table'
import InventoryStockDetailsPanel from '~/components/inventory/stock-details/InventoryStockDetailsPanel.vue'
import InventoryStockTable from '~/components/inventory/InventoryStockTable.vue'
import { useInventoryStockPageQuery } from '~/composables/inventory/useInventoryStockPageQuery'
import { useInventoryStockQuery } from '~/composables/inventory/useInventoryStockQuery'
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

const stockTablePreferenceStore = useInventoryStockTablePreferenceStore()

const selectedStockId = ref<number | null>(props.initialStockId ?? null)
const paginationState = ref<PaginationState>({
  pageIndex: 0,
  pageSize: 50,
})

const effectivePreset = computed<InventoryStockPreset>(() => {
  return props.preset ?? stockTablePreferenceStore.presetState
})
const stockPageQueryParams = computed(() => {
  return {
    preset: effectivePreset.value,
    page: paginationState.value.pageIndex + 1,
    pageSize: paginationState.value.pageSize,
    search: stockTablePreferenceStore.globalFilterState,
    sorting: stockTablePreferenceStore.sortingState,
  }
})
const stocksQuery = useInventoryStockPageQuery(stockPageQueryParams)
const selectedStockQueryId = computed<number>(() => selectedStockId.value ?? 0)
const selectedStockQuery = useInventoryStockQuery(selectedStockQueryId)

const stocks = computed<InventoryStockListItem[]>(() => {
  return stocksQuery.data.value?.results ?? []
})

const totalRowCount = computed<number>(() => {
  return stocksQuery.data.value?.count ?? 0
})

const selectedStock = computed<InventoryStockListItem | null>(() => {
  if (!selectedStockId.value) {
    return null
  }

  for (const stock of stocks.value) {
    if (stock.id === selectedStockId.value) {
      return stock
    }
  }

  return selectedStockQuery.data.value ?? null
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

const setSelectedStockId = (stockId: number | null): void => {
  selectedStockId.value = stockId
  emit('update:selected-stock-id', stockId)
}

/**
 * Opens the stock detail drawer for one selected row.
 *
 * Input example:
 * - `{ id: 14, material: { id: 3, product_name: 'PBS buffer' }, ... }`
 */
const openStockDetails = (stock: InventoryStockListItem): void => {
  setSelectedStockId(stock.id)
}

const closeStockDetails = (): void => {
  setSelectedStockId(null)
}

const updatePaginationState = (nextPaginationState: PaginationState): void => {
  paginationState.value = nextPaginationState
}

const updateSortingState = async (nextSortingState: SortingState): Promise<void> => {
  paginationState.value = {
    ...paginationState.value,
    pageIndex: 0,
  }
  await stockTablePreferenceStore.updateSortingState(nextSortingState)
}

const updateGlobalFilterState = async (nextGlobalFilterState: string): Promise<void> => {
  paginationState.value = {
    ...paginationState.value,
    pageIndex: 0,
  }
  await stockTablePreferenceStore.updateGlobalFilterState(nextGlobalFilterState)
}

watch(
  () => props.preset,
  (nextPreset) => {
    if (!nextPreset) {
      return
    }

    paginationState.value = {
      ...paginationState.value,
      pageIndex: 0,
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
        <p v-else-if="stocks.length === 0" class="text-sm text-slate-600">
          {{ t('inventory.stock_workspace.empty') }}
        </p>
        <InventoryStockTable
          v-else
          :stocks="stocks"
          :pagination-state="paginationState"
          :total-row-count="totalRowCount"
          :sorting-state="stockTablePreferenceStore.sortingState"
          :global-filter-state="stockTablePreferenceStore.globalFilterState"
          :column-order-state="stockTablePreferenceStore.columnOrderState"
          :column-visibility-state="stockTablePreferenceStore.columnVisibilityState"
          @pagination-change="updatePaginationState"
          @sorting-change="updateSortingState"
          @global-filter-change="updateGlobalFilterState"
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
