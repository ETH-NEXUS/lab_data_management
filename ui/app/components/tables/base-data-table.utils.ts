import type { Column, FilterFn, Row } from '@tanstack/vue-table'

export type TableRow = Record<string, unknown>

export type DataTableSemanticGroup = 'operational' | 'identity' | 'lifecycle'

export type DataTableColumnMeta = {
  semanticGroup?: DataTableSemanticGroup
}

/**
 * Applies a TanStack updater function or direct value to state.
 *
 * Accepted data examples:
 * - `updater = [{ id: 'name', desc: false }]`
 * - `updater = (oldState) => [...oldState, nextItem]`
 *
 * Returned data example:
 * - updated state value of type `T`
 */
export const applyUpdater = <T>(updater: T | ((oldState: T) => T), previous: T): T => {
  if (typeof updater === 'function') {
    return (updater as (oldState: T) => T)(previous)
  }
  return updater
}

/**
 * Multi-select column filter that supports scalar and array cell values.
 *
 * Input examples:
 * - scalar row value: `'In stock'`
 * - array row value: `['PCR', 'Cold storage']`
 * - selected filter values: `['in stock']`
 */
export const multiValueFilter: FilterFn<TableRow> = (row, columnId, filterValue) => {
  if (!Array.isArray(filterValue) || filterValue.length === 0) {
    return true
  }

  const selected = new Set(filterValue.map((value) => String(value).toLowerCase()))
  const rawValue = row.getValue(columnId)

  if (Array.isArray(rawValue)) {
    return rawValue.some((entry) => selected.has(String(entry).toLowerCase()))
  }

  return selected.has(String(rawValue ?? '').toLowerCase())
}

/**
 * Collects distinct, non-empty string values for a column filter menu.
 *
 * Returned data example:
 * - `['Cold room A', 'Freezer 2', 'Shelf B3']`
 */
export const getUniqueColumnValues = (rows: Row<TableRow>[], column: Column<TableRow, unknown>): string[] => {
  const values = new Set<string>()

  for (const row of rows) {
    const rawValue = row.getValue(column.id)

    if (Array.isArray(rawValue)) {
      for (const value of rawValue) {
        if (value != null) {
          const token = String(value).trim()
          if (token !== '') {
            values.add(token)
          }
        }
      }
      continue
    }

    if (rawValue != null) {
      const token = String(rawValue).trim()
      if (token !== '') {
        values.add(token)
      }
    }
  }

  return Array.from(values).sort((left, right) => left.localeCompare(right))
}

/**
 * Applies case-insensitive search over filter options.
 *
 * Input example:
 * - `options = ['In stock', 'Low stock']`
 * - `query = 'low'`
 *
 * Returned data example:
 * - `['Low stock']`
 */
export const filterVisibleOptions = (options: string[], query: string): string[] => {
  const normalizedQuery = query.trim().toLowerCase()
  if (normalizedQuery === '') {
    return options
  }
  return options.filter((option) => option.toLowerCase().includes(normalizedQuery))
}

/**
 * Converts a raw cell value to user-facing string for default cell rendering.
 *
 * Input examples:
 * - `['PCR', 'RNA']`
 * - `null`
 * - `12.5`
 *
 * Returned examples:
 * - `'PCR, RNA'`
 * - `'—'`
 * - `'12.5'`
 */
export const formatCellValue = (value: unknown): string => {
  if (Array.isArray(value)) {
    return value.map((item) => String(item)).join(', ')
  }
  if (value == null) {
    return '—'
  }
  return String(value)
}

/**
 * Reads semantic group metadata from a column definition.
 *
 * Accepted metadata example:
 * - `{ semanticGroup: 'identity' }`
 *
 * Returned values:
 * - `'identity' | 'lifecycle' | 'operational' | undefined`
 */
export const getColumnSemanticGroup = (column: Column<TableRow, unknown>): DataTableSemanticGroup | undefined => {
  const columnMeta = column.columnDef.meta as DataTableColumnMeta | undefined
  return columnMeta?.semanticGroup
}
