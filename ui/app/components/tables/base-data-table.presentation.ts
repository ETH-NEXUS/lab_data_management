import type { Column } from '@tanstack/vue-table'
import { getColumnSemanticGroup, type TableRow } from '~/components/tables/base-data-table.utils'

/**
 * Returns the visible label for one column in the column visibility menu.
 *
 * Input example:
 * - `{ id: 'product_name', columnDef: { header: 'Product Name' } }`
 */
export const getColumnVisibilityLabel = (column: Column<TableRow, unknown>): string => {
  const headerValue = column.columnDef.header

  if (typeof headerValue === 'string') {
    return headerValue
  }

  return column.id
}

/**
 * Returns the matching sort icon name for a given sorting state.
 *
 * Input examples:
 * - `false`
 * - `'asc'`
 * - `'desc'`
 */
export const getSortIcon = (state: false | 'asc' | 'desc'): string => {
  if (state === 'asc') return 'i-heroicons-bars-arrow-up'
  if (state === 'desc') return 'i-heroicons-bars-arrow-down'
  return 'i-heroicons-arrows-up-down'
}

/**
 * Maps semantic group metadata to cell background helper classes.
 *
 * Returned value examples:
 * - `'worksheet-cell-group-identity'`
 * - `'worksheet-cell-group-lifecycle'`
 */
export const getSemanticCellClass = (column: Column<TableRow, unknown>): string => {
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
export const getSemanticHeaderClass = (column: Column<TableRow, unknown>): string => {
  const semanticGroup = getColumnSemanticGroup(column)

  if (semanticGroup === 'identity') {
    return 'worksheet-header-group-identity'
  }

  if (semanticGroup === 'lifecycle') {
    return 'worksheet-header-group-lifecycle'
  }

  return 'bg-slate-50'
}

export const getCellClass = (column: Column<TableRow, unknown>): string[] => {
  const classes = ['worksheet-cell', 'px-3 py-2.5 text-sm']
  classes.push(getSemanticCellClass(column))

  return classes
}

export const getHeaderCellClass = (column: Column<TableRow, unknown>): string[] => {
  const classes = ['worksheet-header-cell', 'px-3 py-2.5 text-sm', 'font-semibold text-slate-800']
  classes.push(getSemanticHeaderClass(column))

  return classes
}

export const hasCustomCellRenderer = (column: Column<TableRow, unknown>): boolean => {
  return typeof column.columnDef.cell !== 'undefined'
}

/**
 * Normalizes one optional row class value into a clean string array.
 *
 * Input examples:
 * - `''`
 * - `'text-red-500'`
 * - `['text-red-500', ' font-bold ']`
 *
 * Returned value examples:
 * - `[]`
 * - `['text-red-500']`
 * - `['text-red-500', ' font-bold ']`
 */
export const getRowCellClass = (customClass?: string | string[]): string[] => {
  if (!customClass) {
    return []
  }

  if (Array.isArray(customClass)) {
    return customClass.filter((value) => value.trim() !== '')
  }

  return customClass.trim() === '' ? [] : [customClass]
}
