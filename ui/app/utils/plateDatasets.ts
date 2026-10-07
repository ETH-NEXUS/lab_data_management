import type { Plate } from '~/types/lab'
import { LOG10_SUFFIX } from '~/types/backgroundCorrection'

/**
 * The compact statistics of one measurement of a plate at one time point, as
 * the plate page shows them below the formula of a calculation.
 *
 * Data example:
 * - `{ wells: 64, median: 2.87, min: 1.52, max: 5.1, byWellType: [{ wellType: 'N', median: 3.47 }] }`
 */
export type DatasetSummary = {
  wells: number
  median: number | null
  min: number | null
  max: number | null
  byWellType: { wellType: string; median: number }[]
}

/**
 * The measurement a log10 measurement was calculated from, or null.
 *
 * Examples: `'Lum1_log10'` -> `'Lum1'` (if the plate has `Lum1`); `'Lum1'` -> `null`
 */
export const getLog10SourceLabel = (plate: Plate, label: string | null): string | null => {
  if (!label || !label.endsWith(LOG10_SUFFIX)) return null
  const source = label.slice(0, -LOG10_SUFFIX.length)
  return (plate.details.measurement_labels ?? []).includes(source) ? source : null
}

/**
 * Wells with a value of the source measurement but none of its log10: their
 * value is 0 or below (the server leaves such a well empty in every read).
 */
export const countWellsWithoutLog10 = (plate: Plate, sourceLabel: string, log10Label: string): number => {
  let count = 0
  for (const well of plate.wells ?? []) {
    const sourceValues = well.measurements?.[sourceLabel] ?? []
    const logValues = well.measurements?.[log10Label] ?? []
    if (sourceValues.length > 0 && logValues.length === 0) {
      count += 1
    }
  }
  return count
}

const valueAt = (series: number[] | undefined, index: number): number | null => {
  const value = series?.[index]
  return typeof value === 'number' ? value : null
}

/**
 * Counts, median, min and max of a measurement at one time point, from the
 * statistics the server keeps for every plate.
 */
export const summarizeDataset = (plate: Plate, label: string, timestampIndex: number): DatasetSummary => {
  let wells = 0
  for (const well of plate.wells ?? []) {
    if (valueAt(well.measurements?.[label], timestampIndex) !== null) {
      wells += 1
    }
  }

  const overall = plate.details.overall_stats[label]
  const byWellType: DatasetSummary['byWellType'] = []
  const statsByType = plate.details.stats[label] ?? {}
  for (const wellType of Object.keys(statsByType).sort()) {
    const median = valueAt(statsByType[wellType]?.median, timestampIndex)
    if (median !== null) {
      byWellType.push({ wellType, median })
    }
  }

  return {
    wells,
    median: valueAt(overall?.median, timestampIndex),
    min: valueAt(overall?.min, timestampIndex),
    max: valueAt(overall?.max, timestampIndex),
    byWellType,
  }
}

/**
 * A number for the summary: at most 3 decimals, e.g. `2.86634` -> `'2.866'`, `734` -> `'734'`.
 */
export const formatSummaryNumber = (value: number | null): string => {
  if (value === null) return '–'
  return String(Number(value.toFixed(3)))
}
