import type { PlateDimension } from '~/types/lab'

/**
 * Returns a safe plate identifier string from a dynamic route param.
 *
 * Accepted input examples:
 * - `'42'`
 * - `42`
 * - `['42']`
 * - `undefined`
 *
 * Returned output examples:
 * - `'42'`
 * - `'-'` (fallback when value is missing)
 */
export const getPlateRouteIdLabel = (routeParam: unknown): string => {
  if (Array.isArray(routeParam)) {
    const first = routeParam[0]
    return typeof first === 'string' || typeof first === 'number' ? String(first) : '-'
  }

  if (typeof routeParam === 'string' || typeof routeParam === 'number') {
    return String(routeParam)
  }

  return '-'
}

/**
 * Formats plate barcode for human-readable page titles.
 *
 * Accepted input examples:
 * - `'demo_1'`
 * - `'__TEMPL__library_project_plate'`
 * - `null`
 *
 * Returned output examples:
 * - `'demo_1'`
 * - `'Templates/library/project/plate'`
 * - `'N/A'`
 */
export const formatPlateBarcodeLabel = (barcode: string | null | undefined): string => {
  if (!barcode) return 'N/A'
  if (barcode.startsWith('__TEMPL__')) {
    return `Templates/${barcode.replace('__TEMPL__', '').replaceAll('_', '/')}`
  }
  return barcode
}

/**
 * Converts row/column indices into a numeric well position.
 *
 * Accepted input examples:
 * - `row = 0, col = 0, cols = 24` -> `0`
 * - `row = 1, col = 3, cols = 24` -> `27`
 */
export const positionFromRowCol = (row: number, col: number, dimension: PlateDimension): number => {
  return row * dimension.cols + col
}

/**
 * Converts numeric well position into human-readable plate coordinate (e.g. A1).
 *
 * Accepted input examples:
 * - `position = 0, cols = 24` -> `A1`
 * - `position = 27, cols = 24` -> `B4`
 */
export const hrPositionFromPosition = (position: number, dimension: PlateDimension): string => {
  const row = Math.floor(position / dimension.cols)
  const col = position - row * dimension.cols + 1
  return `${intToChar(row)}${col}`
}

const intToChar = (value: number): string => {
  const code = 'A'.charCodeAt(0)
  return String.fromCharCode(code + value)
}
