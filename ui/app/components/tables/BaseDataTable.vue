<script setup lang="ts">
import {
  FlexRender,
  type ColumnDef,
  type ColumnFiltersState,
  type ColumnOrderState,
  type PaginationState,
  type SortingState,
  type VisibilityState,
} from '@tanstack/vue-table'
import {
  formatCellValue,
  type TableRow,
} from '~/components/tables/base-data-table.utils'
import {
  getCellClass as getPresentationCellClass,
  getHeaderCellClass as getPresentationHeaderCellClass,
  getRowCellClass as getPresentationRowCellClass,
  getSortIcon as getPresentationSortIcon,
  hasCustomCellRenderer as hasPresentationCustomCellRenderer,
} from '~/components/tables/base-data-table.presentation'
import BaseDataTableHeaderFilter from '~/components/tables/BaseDataTableHeaderFilter.vue'
import BaseDataTableToolbar from '~/components/tables/BaseDataTableToolbar.vue'
import { useBaseDataTableState } from '~/components/tables/useBaseDataTableState'

type BaseDataTableProps = {
  data: TableRow[]
  columns: ColumnDef<TableRow, unknown>[]
  globalFilterPlaceholder?: string
  rowCellClassName?: (row: TableRow) => string | string[]
  rowClickable?: boolean
  hideToolbar?: boolean
  enablePagination?: boolean
  pageSize?: number
  paginationState?: PaginationState
  totalRowCount?: number
  manualPagination?: boolean
  manualSorting?: boolean
  manualFiltering?: boolean
  enableColumnFilters?: boolean
  sortingState?: SortingState
  globalFilterState?: string
  columnFiltersState?: ColumnFiltersState
  columnOrderState?: ColumnOrderState
  columnVisibilityState?: VisibilityState
}

const props = withDefaults(defineProps<BaseDataTableProps>(), {
  globalFilterPlaceholder: '',
  rowCellClassName: () => '',
  rowClickable: false,
  hideToolbar: false,
  enablePagination: false,
  pageSize: 50,
  paginationState: undefined,
  totalRowCount: undefined,
  manualPagination: false,
  manualSorting: false,
  manualFiltering: false,
  enableColumnFilters: true,
  sortingState: () => [],
  globalFilterState: '',
  columnFiltersState: () => [],
  columnOrderState: () => [],
  columnVisibilityState: () => ({}),
})

const emit = defineEmits<{
  rowClick: [row: TableRow]
  paginationChange: [pagination: PaginationState]
  sortingChange: [sorting: SortingState]
  globalFilterChange: [value: string]
  columnFiltersChange: [filters: ColumnFiltersState]
  columnOrderChange: [columnOrder: ColumnOrderState]
  columnVisibilityChange: [columnVisibility: VisibilityState]
}>()

const { t } = useI18n()
const {
  globalFilter,
  filterSearch,
  table,
  getVisibleFilterOptions,
  initializeColumnFilterDraft,
  isDraftValueSelected,
  setDraftValueSelection,
  applyColumnFilter,
  clearColumnFilter,
  clearAllFilters,
  totalRowsCount,
} = useBaseDataTableState(props, {
  onPaginationChange: (pagination) => emit('paginationChange', pagination),
  onSortingChange: (sorting) => emit('sortingChange', sorting),
  onGlobalFilterChange: (value) => emit('globalFilterChange', value),
  onColumnFiltersChange: (filters) => emit('columnFiltersChange', filters),
  onColumnOrderChange: (columnOrder) => emit('columnOrderChange', columnOrder),
  onColumnVisibilityChange: (columnVisibility) => emit('columnVisibilityChange', columnVisibility),
})

/**
 * Emits one clicked row to parent components.
 *
 * Input example:
 * - `{ id: 14, product_name: 'PBS Buffer', ... }`
 */
const emitRowClick = (row: TableRow): void => {
  emit('rowClick', row)
}

const getSortIcon = getPresentationSortIcon
const getCellClass = getPresentationCellClass
const getHeaderCellClass = getPresentationHeaderCellClass
const hasCustomCellRenderer = hasPresentationCustomCellRenderer

const getRowCellClass = (row: TableRow): string[] => {
  return getPresentationRowCellClass(props.rowCellClassName?.(row))
}
</script>

<template>
  <div class="space-y-4">
    <BaseDataTableToolbar
      v-if="!props.hideToolbar"
      :table="table"
      :global-filter="globalFilter"
      :global-filter-placeholder="props.globalFilterPlaceholder || t('table.general.search_placeholder')"
      :enable-pagination="props.enablePagination"
      :total-rows-count="totalRowsCount"
      @update:global-filter="(value) => (globalFilter = value)"
      @clear-all="clearAllFilters"
    />

    <div class="worksheet-table-scroll overflow-x-auto rounded-lg border border-[var(--app-border)] bg-white shadow-sm">
      <table class="worksheet-table-grid min-w-full text-sm">
        <thead class="border-b border-[var(--app-border)]">
          <tr v-for="headerGroup in table.getHeaderGroups()" :key="headerGroup.id">
            <th
              v-for="header in headerGroup.headers"
              :key="header.id"
              scope="col"
              :class="getHeaderCellClass(header.column)"
            >
              <template v-if="!header.isPlaceholder">
                <div class="flex items-center justify-between gap-2">
                  <button
                    v-if="header.column.getCanSort()"
                    type="button"
                    class="inline-flex items-center gap-1 rounded-md px-1 py-0.5 text-left text-slate-800 hover:bg-slate-100"
                    @click="header.column.toggleSorting(header.column.getIsSorted() === 'asc')"
                  >
                    <FlexRender :render="header.column.columnDef.header" :props="header.getContext()" />
                    <UIcon :name="getSortIcon(header.column.getIsSorted())" class="h-4 w-4 text-slate-500" />
                  </button>
                  <div v-else class="inline-flex items-center px-1 py-0.5">
                    <FlexRender :render="header.column.columnDef.header" :props="header.getContext()" />
                  </div>

                  <BaseDataTableHeaderFilter
                    v-if="props.enableColumnFilters"
                    :column="header.column"
                    :filter-search="filterSearch"
                    :get-visible-filter-options="getVisibleFilterOptions"
                    :initialize-column-filter-draft="initializeColumnFilterDraft"
                    :is-draft-value-selected="isDraftValueSelected"
                    :set-draft-value-selection="setDraftValueSelection"
                    :apply-column-filter="applyColumnFilter"
                    :clear-column-filter="clearColumnFilter"
                  />
                </div>
              </template>
            </th>
          </tr>
        </thead>

        <tbody>
          <tr
            v-for="row in table.getRowModel().rows"
            :key="row.id"
            :class="[
              'worksheet-row border-b border-slate-100 transition-colors hover:bg-slate-50',
              props.rowClickable ? 'cursor-pointer' : '',
            ]"
            @click="props.rowClickable ? emitRowClick(row.original) : undefined"
          >
            <td
              v-for="cell in row.getVisibleCells()"
              :key="cell.id"
              :class="[...getCellClass(cell.column), 'text-slate-700', ...getRowCellClass(row.original)]"
            >
              <template v-if="hasCustomCellRenderer(cell.column)">
                <FlexRender :render="cell.column.columnDef.cell" :props="cell.getContext()" />
              </template>
              <template v-else>
                {{ formatCellValue(cell.getValue()) }}
              </template>
            </td>
          </tr>

          <tr v-if="table.getRowModel().rows.length === 0">
            <td :colspan="table.getVisibleLeafColumns().length" class="px-3 py-10 text-center text-sm text-slate-500">
              {{ t('table.general.no_rows') }}
            </td>
          </tr>
        </tbody>
      </table>
    </div>
  </div>
</template>
