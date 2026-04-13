<script setup lang="ts">
import {
  FlexRender,
  getCoreRowModel,
  getFilteredRowModel,
  getPaginationRowModel,
  getSortedRowModel,
  useVueTable,
  type Column,
  type ColumnDef,
  type ColumnFiltersState,
  type PaginationState,
  type SortingState,
} from '@tanstack/vue-table'
import {
  applyUpdater,
  buildFrozenLeftByColumnId,
  filterVisibleOptions,
  formatCellValue,
  getColumnSemanticGroup,
  getUniqueColumnValues,
  multiValueFilter,
  type TableRow,
} from '~/components/tables/base-data-table.utils'

type BaseDataTableProps = {
  data: TableRow[]
  columns: ColumnDef<TableRow, unknown>[]
  globalFilterPlaceholder?: string
  frozenColumnCount?: number
  rowCellClassName?: (row: TableRow) => string | string[]
  rowClickable?: boolean
  hideToolbar?: boolean
  enablePagination?: boolean
  pageSize?: number
}

const props = withDefaults(defineProps<BaseDataTableProps>(), {
  globalFilterPlaceholder: '',
  frozenColumnCount: 0,
  rowCellClassName: () => '',
  rowClickable: false,
  hideToolbar: false,
  enablePagination: false,
  pageSize: 50,
})

const emit = defineEmits<{
  rowClick: [row: TableRow]
}>()

const { t } = useI18n()

const sorting = ref<SortingState>([])
const globalFilter = ref('')
const columnFilters = ref<ColumnFiltersState>([])
const pagination = ref<PaginationState>({
  pageIndex: 0,
  pageSize: props.pageSize,
})

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
    get pagination() {
      return pagination.value
    },
  },
  filterFns: {
    multiValue: multiValueFilter,
  },
  onSortingChange: (updater) => {
    sorting.value = applyUpdater(updater, sorting.value)
  },
  onGlobalFilterChange: (updater) => {
    globalFilter.value = applyUpdater(updater, globalFilter.value)
  },
  onColumnFiltersChange: (updater) => {
    columnFilters.value = applyUpdater(updater, columnFilters.value)
  },
  onPaginationChange: (updater) => {
    pagination.value = applyUpdater(updater, pagination.value)
  },
  getCoreRowModel: getCoreRowModel(),
  getFilteredRowModel: getFilteredRowModel(),
  getSortedRowModel: getSortedRowModel(),
  getPaginationRowModel: props.enablePagination ? getPaginationRowModel() : undefined,
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

/**
 * Filters the option list shown in one column filter popover.
 *
 * Input example:
 * - `column.id = 'item_type'`
 *
 * Returned example:
 * - `['Consumable', 'Reagent']`
 */
const getVisibleFilterOptions = (column: Column<TableRow, unknown>): string[] => {
  return filterVisibleOptions(uniqueColumnValues.value[column.id] ?? [], filterSearch.value[column.id] ?? '')
}

/**
 * Initializes temporary draft values when one filter popover opens.
 *
 * Input example:
 * - `column.id = 'inventory_status'`
 */
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

/**
 * Applies one column filter using the selected draft values.
 *
 * Input example:
 * - `column.id = 'location'`
 * - `filterDraft['location'] = ['Room A / Shelf 2']`
 */
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

/**
 * Returns the matching sort icon name for a given sorting state.
 *
 * Input examples:
 * - `false`
 * - `'asc'`
 * - `'desc'`
 */
const getSortIcon = (state: false | 'asc' | 'desc'): string => {
  if (state === 'asc') return 'i-heroicons-bars-arrow-up'
  if (state === 'desc') return 'i-heroicons-bars-arrow-down'
  return 'i-heroicons-arrows-up-down'
}

const frozenLeafColumns = computed(() => {
  return table.getAllLeafColumns().slice(0, props.frozenColumnCount)
})

const frozenLeftByColumnId = computed<Record<string, number>>(() => {
  return buildFrozenLeftByColumnId(frozenLeafColumns.value)
})

const isFrozenColumn = (column: Column<TableRow, unknown>): boolean => {
  return props.frozenColumnCount > 0 && frozenLeftByColumnId.value[column.id] != null
}

/**
 * Maps semantic group metadata to cell background helper classes.
 *
 * Returned value examples:
 * - `'worksheet-cell-group-identity'`
 * - `'worksheet-cell-group-lifecycle'`
 */
const getSemanticCellClass = (column: Column<TableRow, unknown>): string => {
  const semanticGroup = getColumnSemanticGroup(column)

  if (semanticGroup === 'identity') {
    return 'worksheet-cell-group-identity'
  }

  if (semanticGroup === 'lifecycle') {
    return 'worksheet-cell-group-lifecycle'
  }

  return 'bg-white'
}

/**
 * Maps semantic group metadata to header background helper classes.
 *
 * Returned value examples:
 * - `'worksheet-header-group-identity'`
 * - `'worksheet-header-group-lifecycle'`
 */
const getSemanticHeaderClass = (column: Column<TableRow, unknown>): string => {
  const semanticGroup = getColumnSemanticGroup(column)

  if (semanticGroup === 'identity') {
    return 'worksheet-header-group-identity'
  }

  if (semanticGroup === 'lifecycle') {
    return 'worksheet-header-group-lifecycle'
  }

  return 'bg-slate-50'
}

const getCellClass = (column: Column<TableRow, unknown>): string[] => {
  const classes = ['worksheet-cell', 'px-3 py-2.5 text-sm']

  if (isFrozenColumn(column)) {
    // Frozen columns always use white background for scan stability.
    classes.push('sticky z-20 bg-white worksheet-sticky-cell')
  } else {
    classes.push(getSemanticCellClass(column))
  }

  return classes
}

const getHeaderCellClass = (column: Column<TableRow, unknown>): string[] => {
  const classes = ['worksheet-header-cell', 'px-3 py-2.5 text-sm', 'font-semibold text-slate-800']

  if (isFrozenColumn(column)) {
    // Frozen columns always use white background for scan stability.
    classes.push('sticky z-30 bg-white worksheet-sticky-header')
  } else {
    classes.push(getSemanticHeaderClass(column))
  }

  return classes
}

/**
 * Returns sticky position styles for frozen columns.
 *
 * Returned style example:
 * - `{ left: '184px', minWidth: '220px', width: '220px' }`
 */
const getStickyStyle = (column: Column<TableRow, unknown>): Record<string, string> => {
  if (!isFrozenColumn(column)) {
    return {}
  }

  return {
    left: `${frozenLeftByColumnId.value[column.id]}px`,
    minWidth: `${column.getSize()}px`,
    width: `${column.getSize()}px`,
    boxShadow: '1px 0 0 rgba(148, 163, 184, 0.35)',
  }
}

/**
 * Emits one clicked row to parent components.
 *
 * Input example:
 * - `{ id: 14, product_name: 'PBS Buffer', ... }`
 */
const emitRowClick = (row: TableRow): void => {
  emit('rowClick', row)
}

const hasCustomCellRenderer = (column: Column<TableRow, unknown>): boolean => {
  return typeof column.columnDef.cell !== 'undefined'
}

const getRowCellClass = (row: TableRow): string[] => {
  const customClass = props.rowCellClassName?.(row)

  if (!customClass) {
    return []
  }

  if (Array.isArray(customClass)) {
    return customClass.filter((value) => value.trim() !== '')
  }

  return customClass.trim() === '' ? [] : [customClass]
}

const totalRowsCount = computed<number>(() => {
  return props.enablePagination ? table.getFilteredRowModel().rows.length : table.getRowModel().rows.length
})
</script>

<template>
  <div class="space-y-4">
    <div v-if="!props.hideToolbar" class="flex flex-col gap-2 sm:flex-row sm:items-center sm:justify-between">
      <UInput
        v-model="globalFilter"
        :placeholder="props.globalFilterPlaceholder || t('table.general.search_placeholder')"
        icon="i-heroicons-magnifying-glass"
        class="w-full sm:max-w-sm"
      />

      <div class="flex items-center gap-2">
        <UBadge color="neutral" variant="soft" class="font-medium">
          {{ t('table.general.rows', { count: totalRowsCount }) }}
        </UBadge>
        <div v-if="props.enablePagination && totalRowsCount > 0" class="flex items-center gap-1">
          <span class="text-xs text-slate-600">
            {{ table.getState().pagination.pageIndex + 1 }} / {{ table.getPageCount() }}
          </span>
          <UButton
            size="xs"
            color="neutral"
            variant="ghost"
            icon="i-heroicons-chevron-left"
            :disabled="!table.getCanPreviousPage()"
            @click="table.previousPage()"
          />
          <UButton
            size="xs"
            color="neutral"
            variant="ghost"
            icon="i-heroicons-chevron-right"
            :disabled="!table.getCanNextPage()"
            @click="table.nextPage()"
          />
        </div>
        <UButton
          variant="ghost"
          color="neutral"
          icon="i-heroicons-x-mark"
          :label="t('table.general.clear_all')"
          @click="clearAllFilters"
        />
      </div>
    </div>

    <div class="worksheet-table-scroll overflow-x-auto rounded-lg border border-[var(--app-border)] bg-white shadow-sm">
      <table class="worksheet-table-grid min-w-full text-sm">
        <thead class="border-b border-[var(--app-border)]">
          <tr v-for="headerGroup in table.getHeaderGroups()" :key="headerGroup.id">
            <th
              v-for="header in headerGroup.headers"
              :key="header.id"
              scope="col"
              :class="getHeaderCellClass(header.column)"
              :style="getStickyStyle(header.column)"
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

                  <UPopover v-if="header.column.getCanFilter()" :content="{ side: 'bottom', align: 'end' }">
                    <UButton
                      variant="ghost"
                      color="neutral"
                      size="xs"
                      icon="i-heroicons-funnel"
                      :class="header.column.getFilterValue() ? 'text-blue-700' : 'text-slate-500'"
                      @click="initializeColumnFilterDraft(header.column)"
                    />

                    <template #content="{ close }">
                      <div class="w-64 space-y-3 p-3">
                        <div class="flex flex-col gap-1">
                          <UButton
                            variant="ghost"
                            color="neutral"
                            icon="i-heroicons-bars-arrow-up"
                            :label="t('table.general.sort_asc')"
                            block
                            @click="
                              () => {
                                header.column.toggleSorting(false)
                                close()
                              }
                            "
                          />
                          <UButton
                            variant="ghost"
                            color="neutral"
                            icon="i-heroicons-bars-arrow-down"
                            :label="t('table.general.sort_desc')"
                            block
                            @click="
                              () => {
                                header.column.toggleSorting(true)
                                close()
                              }
                            "
                          />
                          <UButton
                            variant="ghost"
                            color="neutral"
                            icon="i-heroicons-arrows-up-down"
                            :label="t('table.general.clear_sort')"
                            block
                            @click="
                              () => {
                                header.column.clearSorting()
                                close()
                              }
                            "
                          />
                        </div>

                        <UDivider />

                        <UInput
                          v-model="filterSearch[header.column.id]"
                          :placeholder="t('table.general.filter_options_placeholder')"
                          icon="i-heroicons-magnifying-glass"
                          size="sm"
                        />

                        <div class="max-h-44 space-y-1 overflow-y-auto pr-1">
                          <label
                            v-for="option in getVisibleFilterOptions(header.column)"
                            :key="`${header.column.id}-${option}`"
                            class="flex cursor-pointer items-center gap-2 rounded px-1 py-1 hover:bg-slate-50"
                          >
                            <UCheckbox
                              :model-value="isDraftValueSelected(header.column.id, option)"
                              @update:model-value="
                                (checked) => setDraftValueSelection(header.column.id, option, Boolean(checked))
                              "
                            />
                            <span class="truncate text-xs text-slate-700">{{ option }}</span>
                          </label>

                          <p
                            v-if="getVisibleFilterOptions(header.column).length === 0"
                            class="px-1 py-1 text-xs text-slate-500"
                          >
                            {{ t('table.general.no_filter_options') }}
                          </p>
                        </div>

                        <div class="flex items-center gap-2">
                          <UButton
                            color="primary"
                            size="sm"
                            :label="t('table.general.apply_filter')"
                            @click="applyColumnFilter(header.column, close)"
                          />
                          <UButton
                            variant="ghost"
                            color="neutral"
                            size="sm"
                            :label="t('table.general.clear_filter')"
                            @click="clearColumnFilter(header.column, close)"
                          />
                        </div>
                      </div>
                    </template>
                  </UPopover>
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
              v-for="cell in row.getAllCells()"
              :key="cell.id"
              :class="[...getCellClass(cell.column), 'text-slate-700', ...getRowCellClass(row.original)]"
              :style="getStickyStyle(cell.column)"
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
            <td :colspan="table.getAllLeafColumns().length" class="px-3 py-10 text-center text-sm text-slate-500">
              {{ t('table.general.no_rows') }}
            </td>
          </tr>
        </tbody>
      </table>
    </div>
  </div>
</template>
