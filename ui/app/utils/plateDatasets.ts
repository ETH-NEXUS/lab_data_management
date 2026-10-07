import type { Plate } from '~/types/lab'
import { LOG10_SUFFIX, type ActivityDataset } from '~/types/backgroundCorrection'
import { getStatsSeriesValue } from '~/utils/plateStats'

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

/**
 * Counts, median, min and max of a measurement at one time point, from the
 * statistics the server keeps for every plate.
 */
export const summarizeDataset = (plate: Plate, label: string, timestampIndex: number): DatasetSummary => {
  let wells = 0
  for (const well of plate.wells ?? []) {
    if (getStatsSeriesValue(well.measurements?.[label], timestampIndex) !== null) {
      wells += 1
    }
  }

  const overall = plate.details.overall_stats[label]
  const byWellType: DatasetSummary['byWellType'] = []
  const statsByType = plate.details.stats[label] ?? {}
  for (const wellType of Object.keys(statsByType).sort()) {
    const median = getStatsSeriesValue(statsByType[wellType]?.median, timestampIndex)
    if (median !== null) {
      byWellType.push({ wellType, median })
    }
  }

  return {
    wells,
    median: getStatsSeriesValue(overall?.median, timestampIndex),
    min: getStatsSeriesValue(overall?.min, timestampIndex),
    max: getStatsSeriesValue(overall?.max, timestampIndex),
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

/**
 * How a %Activity measurement was calculated, from its name, or null.
 * Every measurement and pair of its well types could be the source, so they
 * are looked up by the name they give.
 *
 * Example: `'Lum1_activity_N_P'` -> `{ label: 'Lum1_activity_N_P', source: 'Lum1', negativeType: 'N', positiveType: 'P' }`
 */
export const getActivityDataset = (plate: Plate, label: string | null): ActivityDataset | null => {
  if (!label || !label.includes('_activity_')) return null

  for (const source of plate.details.measurement_labels ?? []) {
    const wellTypes = Object.keys(plate.details.stats[source] ?? {})
    for (const negativeType of wellTypes) {
      for (const positiveType of wellTypes) {
        if (label === `${source}_activity_${negativeType}_${positiveType}`) {
          return { label, source, negativeType, positiveType }
        }
      }
    }
  }
  return null
}

/**
 * The median of one well type of a measurement at one time point, e.g. of the N wells of Lum1.
 */
export const getWellTypeMedian = (
  plate: Plate,
  label: string,
  wellType: string,
  timestampIndex: number,
): number | null => {
  return getStatsSeriesValue(plate.details.stats[label]?.[wellType]?.median, timestampIndex)
}
