import type { Plate } from '~/types/lab'
import {
  ACTIVITY_INFIX,
  BACKGROUND_CORRECTION_METHODS,
  CORRECTION_INFIX,
  LOG10_SUFFIX,
  type ActivityDataset,
  type BackgroundCorrection,
} from '~/types/plateCalculations'
import { getStatsSeriesValue } from '~/utils/plateStats'

/**
 * The measurements of a plate, and how a calculated one was calculated: the
 * server names a new measurement after its calculation (see types/plateCalculations.ts),
 * so the plate page finds the calculation by that name.
 */

/**
 * Well types that have a value of the measurement on the plate, e.g. `['C', 'N1', 'P']`.
 */
export const getWellTypesOfMeasurement = (plate: Plate, label: string | null): string[] => {
  if (!label) return []
  return Object.keys(plate.details.stats[label] ?? {}).sort()
}

// The name the server gives a background correction, e.g. ('Lum1', 'R', 'median') -> 'Lum1_bc_R_median'
const correctionLabel = (source: string, referenceType: string, method: string): string =>
  `${source}${CORRECTION_INFIX}${referenceType}_${method}`

/**
 * Background corrections of one measurement of the plate. Each well type of the
 * measurement could have been the reference, so every well type and method is
 * looked up by the name it gives.
 *
 * Returned data example (plate with `Lum1_bc_R_median` and `Lum1_bc_R_mean`):
 * - `[{ label: 'Lum1_bc_R_median', source: 'Lum1', referenceType: 'R', method: 'median' },
 *     { label: 'Lum1_bc_R_mean', source: 'Lum1', referenceType: 'R', method: 'mean' }]`
 */
export const findBackgroundCorrections = (plate: Plate, label: string | null): BackgroundCorrection[] => {
  if (!label) return []

  const plateLabels = plate.details.measurement_labels ?? []
  const corrections: BackgroundCorrection[] = []
  for (const referenceType of getWellTypesOfMeasurement(plate, label)) {
    for (const method of BACKGROUND_CORRECTION_METHODS) {
      const corrected = correctionLabel(label, referenceType, method)
      if (plateLabels.includes(corrected)) {
        corrections.push({ label: corrected, source: label, referenceType, method })
      }
    }
  }
  return corrections
}

/**
 * How a background correction was calculated, or null if the measurement is none.
 *
 * Example: `'Lum1_bc_R_median'` -> `{ label: 'Lum1_bc_R_median', source: 'Lum1', referenceType: 'R', method: 'median' }`
 */
export const getCorrectionDataset = (plate: Plate, label: string | null): BackgroundCorrection | null => {
  if (!label || !label.includes(CORRECTION_INFIX)) return null

  for (const source of plate.details.measurement_labels ?? []) {
    const found = findBackgroundCorrections(plate, source).find((correction) => correction.label === label)
    if (found) return found
  }
  return null
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
 * How a %Activity measurement was calculated, or null. Every measurement and
 * pair of its well types could be the source, so they are looked up by the
 * name they give.
 *
 * Example: `'Lum1_activity_N_P'` -> `{ label: 'Lum1_activity_N_P', source: 'Lum1', negativeType: 'N', positiveType: 'P' }`
 */
export const getActivityDataset = (plate: Plate, label: string | null): ActivityDataset | null => {
  if (!label || !label.includes(ACTIVITY_INFIX)) return null

  for (const source of plate.details.measurement_labels ?? []) {
    const wellTypes = getWellTypesOfMeasurement(plate, source)
    for (const negativeType of wellTypes) {
      for (const positiveType of wellTypes) {
        if (label === `${source}${ACTIVITY_INFIX}${negativeType}_${positiveType}`) {
          return { label, source, negativeType, positiveType }
        }
      }
    }
  }
  return null
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
  for (const wellType of getWellTypesOfMeasurement(plate, label)) {
    const median = getWellTypeStatistic(plate, label, wellType, 'median', timestampIndex)
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
 * The median or mean of one well type of a measurement at one time point,
 * e.g. the median of the R wells of Lum1.
 */
export const getWellTypeStatistic = (
  plate: Plate,
  label: string,
  wellType: string,
  statistic: 'median' | 'mean',
  timestampIndex: number,
): number | null => {
  return getStatsSeriesValue(plate.details.stats[label]?.[wellType]?.[statistic], timestampIndex)
}

/**
 * A number for the summary: at most 3 decimals, e.g. `2.86634` -> `'2.866'`, `734` -> `'734'`.
 */
export const formatSummaryNumber = (value: number | null): string => {
  if (value === null) return '–'
  return String(Number(value.toFixed(3)))
}
