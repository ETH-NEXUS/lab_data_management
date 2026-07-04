<script setup lang="ts">
import { computed } from 'vue'
import type {
  ColumnFiltersState,
  ColumnOrderState,
  PaginationState,
  SortingState,
  VisibilityState,
} from '@tanstack/vue-table'
import BaseDataTable from '~/components/tables/BaseDataTable.vue'
import { createInventoryStockTableColumns } from '~/components/inventory/inventory-stock-table.columns'
import type { InventoryStockListItem } from '~/types/inventory'

type Props = {
  stocks: InventoryStockListItem[]
  paginationState?: PaginationState
  totalRowCount?: number
  sortingState?: SortingState
  globalFilterState?: string
  columnFiltersState?: ColumnFiltersState
  columnOrderState?: ColumnOrderState
  columnVisibilityState?: VisibilityState
}

const props = defineProps<Props>()

const emit = defineEmits<{
  (e: 'select-stock', stock: InventoryStockListItem): void
  (e: 'pagination-change', pagination: PaginationState): void
  (e: 'sorting-change', sorting: SortingState): void
  (e: 'global-filter-change', value: string): void
  (e: 'column-filters-change', filters: ColumnFiltersState): void
  (e: 'column-order-change', columnOrder: ColumnOrderState): void
  (e: 'column-visibility-change', columnVisibility: VisibilityState): void
}>()

const { t } = useI18n()

const onSelectStock = (stock: InventoryStockListItem): void => {
  emit('select-stock', stock)
}

/**
 * Returns custom table cell class for one row.
 *
 * Input example:
 * - `{ id: 14, is_favorite: true, ... }`
 *
 * Returned example:
 * - `'!bg-amber-100'`
 */
const getStockRowCellClass = (row: Record<string, unknown>): string => {
  const stock = row as InventoryStockListItem

  if (stock.is_favorite) {
    return '!bg-amber-100'
  }
  return ''
}

/**
 * Builds columns with a product-name click action that opens stock details.
 *
 * Input callback example:
 * - `(stock) => emit('select-stock', stock)`
 */
const tableColumns = computed(() => createInventoryStockTableColumns(t, onSelectStock))
</script>

<template>
  <BaseDataTable
    :data="props.stocks as unknown as Record<string, unknown>[]"
    :columns="tableColumns"
    :row-cell-class-name="getStockRowCellClass"
    enable-pagination
    :page-size="50"
    :pagination-state="props.paginationState"
    :total-row-count="props.totalRowCount"
    manual-pagination
    :sorting-state="props.sortingState"
    :global-filter-state="props.globalFilterState"
    :column-filters-state="props.columnFiltersState"
    :column-order-state="props.columnOrderState"
    :column-visibility-state="props.columnVisibilityState"
    :global-filter-placeholder="t('table.general.search_placeholder')"
    @pagination-change="(pagination) => emit('pagination-change', pagination)"
    @sorting-change="(sorting) => emit('sorting-change', sorting)"
    @global-filter-change="(value) => emit('global-filter-change', value)"
    @column-filters-change="(filters) => emit('column-filters-change', filters)"
    @column-order-change="(columnOrder) => emit('column-order-change', columnOrder)"
    @column-visibility-change="(columnVisibility) => emit('column-visibility-change', columnVisibility)"
  />
</template>
