import { computed, ref, watch, type ComputedRef, type Ref } from 'vue'
import {
  getCoreRowModel,
  getFilteredRowModel,
  getPaginationRowModel,
  getSortedRowModel,
  useVueTable,
  type Column,
  type ColumnDef,
  type ColumnFiltersState,
  type ColumnOrderState,
  type PaginationState,
  type SortingState,
  type VisibilityState,
} from '@tanstack/vue-table'
import {
  applyUpdater,
  filterVisibleOptions,
  getUniqueColumnValues,
  multiValueFilter,
  type TableRow,
} from '~/components/tables/base-data-table.utils'

export type BaseDataTableStateProps = {
  data: TableRow[]
  columns: ColumnDef<TableRow, unknown>[]
  enablePagination: boolean
  pageSize: number
  paginationState?: PaginationState
  totalRowCount?: number
  manualPagination?: boolean
  manualSorting?: boolean
  manualFiltering?: boolean
  sortingState: SortingState
  globalFilterState: string
  columnFiltersState: ColumnFiltersState
  columnOrderState: ColumnOrderState
  columnVisibilityState: VisibilityState
}

export type BaseDataTableStateEmitters = {
  onPaginationChange: (pagination: PaginationState) => void
  onSortingChange: (sorting: SortingState) => void
  onGlobalFilterChange: (value: string) => void
  onColumnFiltersChange: (filters: ColumnFiltersState) => void
  onColumnOrderChange: (columnOrder: ColumnOrderState) => void
  onColumnVisibilityChange: (columnVisibility: VisibilityState) => void
}

export type BaseDataTableStateResult = {
  sorting: Ref<SortingState>
  globalFilter: Ref<string>
  columnFilters: Ref<ColumnFiltersState>
  columnOrder: Ref<ColumnOrderState>
  columnVisibility: Ref<VisibilityState>
  pagination: Ref<PaginationState>
  filterDraft: Ref<Record<string, string[]>>
  filterSearch: Ref<Record<string, string>>
  table: ReturnType<typeof useVueTable<TableRow>>
  uniqueColumnValues: ComputedRef<Record<string, string[]>>
  getVisibleFilterOptions: (column: Column<TableRow, unknown>) => string[]
  initializeColumnFilterDraft: (column: Column<TableRow, unknown>) => void
  isDraftValueSelected: (columnId: string, option: string) => boolean
  setDraftValueSelection: (columnId: string, option: string, checked: boolean) => void
  applyColumnFilter: (column: Column<TableRow, unknown>, close?: () => void) => void
  clearColumnFilter: (column: Column<TableRow, unknown>, close?: () => void) => void
  clearAllFilters: () => void
  totalRowsCount: ComputedRef<number>
}

/**
 * Owns TanStack table state, filter drafts, and controlled-state synchronization
 * for `BaseDataTable.vue`.
 */
export const useBaseDataTableState = (
  props: BaseDataTableStateProps,
  emitters: BaseDataTableStateEmitters,
): BaseDataTableStateResult => {
  const sorting = ref<SortingState>(props.sortingState)
  const globalFilter = ref(props.globalFilterState)
  const columnFilters = ref<ColumnFiltersState>(props.columnFiltersState)
  const columnOrder = ref<ColumnOrderState>(props.columnOrderState)
  const columnVisibility = ref<VisibilityState>(props.columnVisibilityState)
  const pagination = ref<PaginationState>(
    props.paginationState ?? {
      pageIndex: 0,
      pageSize: props.pageSize,
    },
  )

  const filterDraft = ref<Record<string, string[]>>({})
  const filterSearch = ref<Record<string, string>>({})

  const table = useVueTable({
    get data() {
      return props.data
    },
    get columns() {
      return props.columns
    },
    defaultColumn: {
      filterFn: multiValueFilter,
    },
    state: {
      get sorting() {
        return sorting.value
      },
      get globalFilter() {
        return globalFilter.value
      },
      get columnFilters() {
        return columnFilters.value
      },
      get columnOrder() {
        return columnOrder.value
      },
      get columnVisibility() {
        return columnVisibility.value
      },
      get pagination() {
        return pagination.value
      },
    },
    filterFns: {
      multiValue: multiValueFilter,
    },
    manualPagination: props.manualPagination ?? false,
    manualSorting: props.manualSorting ?? false,
    manualFiltering: props.manualFiltering ?? false,
    rowCount: props.totalRowCount,
    onSortingChange: (updater) => {
      sorting.value = applyUpdater(updater, sorting.value)
      emitters.onSortingChange(sorting.value)
    },
    onGlobalFilterChange: (updater) => {
      globalFilter.value = applyUpdater(updater, globalFilter.value)
      emitters.onGlobalFilterChange(globalFilter.value)
    },
    onColumnFiltersChange: (updater) => {
      columnFilters.value = applyUpdater(updater, columnFilters.value)
      emitters.onColumnFiltersChange(columnFilters.value)
    },
    onColumnOrderChange: (updater) => {
      columnOrder.value = applyUpdater(updater, columnOrder.value)
      emitters.onColumnOrderChange(columnOrder.value)
    },
    onColumnVisibilityChange: (updater) => {
      columnVisibility.value = applyUpdater(updater, columnVisibility.value)
      emitters.onColumnVisibilityChange(columnVisibility.value)
    },
    onPaginationChange: (updater) => {
      pagination.value = applyUpdater(updater, pagination.value)
      emitters.onPaginationChange(pagination.value)
    },
    getCoreRowModel: getCoreRowModel(),
    getFilteredRowModel: props.manualFiltering ? undefined : getFilteredRowModel(),
    getSortedRowModel: props.manualSorting ? undefined : getSortedRowModel(),
    getPaginationRowModel:
      props.enablePagination && !(props.manualPagination ?? false) ? getPaginationRowModel() : undefined,
  })

  /**
   * Builds unique column values from unfiltered rows to populate filter popovers.
   *
   * Returned structure example:
   * - `{ inventory_status: ['In stock', 'Low stock'], device_type: ['Freezer'] }`
   */
  const uniqueColumnValues = computed<Record<string, string[]>>(() => {
    const valuesByColumnId: Record<string, string[]> = {}
    const rows = table.getPreFilteredRowModel().rows

    for (const column of table.getAllLeafColumns()) {
      valuesByColumnId[column.id] = getUniqueColumnValues(rows, column)
    }

    return valuesByColumnId
  })

  const getVisibleFilterOptions = (column: Column<TableRow, unknown>): string[] => {
    return filterVisibleOptions(uniqueColumnValues.value[column.id] ?? [], filterSearch.value[column.id] ?? '')
  }

  const initializeColumnFilterDraft = (column: Column<TableRow, unknown>): void => {
    const filterValue = column.getFilterValue()
    filterDraft.value[column.id] = Array.isArray(filterValue) ? filterValue.map((value) => String(value)) : []
    filterSearch.value[column.id] = ''
  }

  const isDraftValueSelected = (columnId: string, option: string): boolean => {
    return (filterDraft.value[columnId] ?? []).includes(option)
  }

  const setDraftValueSelection = (columnId: string, option: string, checked: boolean): void => {
    const current = filterDraft.value[columnId] ?? []
    const next = checked ? Array.from(new Set([...current, option])) : current.filter((value) => value !== option)
    filterDraft.value[columnId] = next
  }

  const applyColumnFilter = (column: Column<TableRow, unknown>, close?: () => void): void => {
    const nextFilter = filterDraft.value[column.id] ?? []
    column.setFilterValue(nextFilter.length > 0 ? nextFilter : undefined)
    close?.()
  }

  const clearColumnFilter = (column: Column<TableRow, unknown>, close?: () => void): void => {
    filterDraft.value[column.id] = []
    filterSearch.value[column.id] = ''
    column.setFilterValue(undefined)
    close?.()
  }

  const clearAllFilters = (): void => {
    sorting.value = []
    globalFilter.value = ''
    columnFilters.value = []
    filterDraft.value = {}
    filterSearch.value = {}
    pagination.value.pageIndex = 0
    emitters.onSortingChange(sorting.value)
    emitters.onGlobalFilterChange(globalFilter.value)
    emitters.onColumnFiltersChange(columnFilters.value)
    emitters.onPaginationChange(pagination.value)
  }

  watch(
    () => props.pageSize,
    (pageSize) => {
      pagination.value = {
        pageIndex: 0,
        pageSize,
      }
    },
  )

  watch(
    () => props.paginationState,
    (nextPaginationState) => {
      if (!nextPaginationState) {
        return
      }

      pagination.value = nextPaginationState
    },
  )

  watch(
    () => props.sortingState,
    (nextSortingState) => {
      sorting.value = nextSortingState
    },
  )

  watch(
    () => props.globalFilterState,
    (nextGlobalFilterState) => {
      globalFilter.value = nextGlobalFilterState
    },
  )

  watch(
    () => props.columnFiltersState,
    (nextColumnFiltersState) => {
      columnFilters.value = nextColumnFiltersState
    },
  )

  watch(
    () => props.columnOrderState,
    (nextColumnOrderState) => {
      columnOrder.value = nextColumnOrderState
    },
  )

  watch(
    () => props.columnVisibilityState,
    (nextColumnVisibilityState) => {
      columnVisibility.value = nextColumnVisibilityState
    },
  )

  const totalRowsCount = computed<number>(() => {
    if ((props.manualPagination ?? false) && typeof props.totalRowCount === 'number') {
      return props.totalRowCount
    }

    return props.enablePagination ? table.getFilteredRowModel().rows.length : table.getRowModel().rows.length
  })

  return {
    sorting,
    globalFilter,
    columnFilters,
    columnOrder,
    columnVisibility,
    pagination,
    filterDraft,
    filterSearch,
    table,
    uniqueColumnValues,
    getVisibleFilterOptions,
    initializeColumnFilterDraft,
    isDraftValueSelected,
    setDraftValueSelection,
    applyColumnFilter,
    clearColumnFilter,
    clearAllFilters,
    totalRowsCount,
  }
}
