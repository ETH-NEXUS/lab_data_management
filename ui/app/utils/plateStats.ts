import type { Plate, PlateStats } from '~/types/lab'
import type { PlateControlLabelOption, PlateMeasurementStatsMap, PlateTimestampOption } from '~/types/plates'

const getMeasurementStats = (plate: Plate, measurement: string | null): PlateMeasurementStatsMap | null => {
  if (!measurement) return null
  return plate.details.stats[measurement] ?? null
}

const getStatsSeriesValue = (series: number[] | undefined, timestampIndex: number): number | null => {
  if (!series || timestampIndex < 0 || timestampIndex >= series.length) return null
  const value = series[timestampIndex]
  return typeof value === 'number' ? value : null
}

const getPlateStatsForLabel = (plate: Plate, measurement: string | null, label: string | null): PlateStats | null => {
  if (!measurement || !label) return null
  const measurementStats = getMeasurementStats(plate, measurement)
  if (!measurementStats) return null
  return measurementStats[label] ?? null
}

/**
 * Builds timestamp options for selected measurement.
 *
 * Accepted data example:
 * - `measurement = 'OD600'`
 *
 * Returned data example:
 * - `[{ label: '2026-01-01T00:00:00Z', value: 0 }]`
 */
export const getTimestampOptions = (plate: Plate, measurement: string | null): PlateTimestampOption[] => {
  if (!measurement) return []
  const timestamps = plate.details.measurement_timestamps[measurement] ?? []
  return timestamps.map((timestamp, index) => ({
    label: timestamp,
    value: index,
  }))
}

/**
 * Builds positive/negative control options from available stats labels.
 *
 * Accepted data example:
 * - stats labels: `['P1', 'N1', 'C']`
 *
 * Returned data example:
 * - `[{ label: 'P1', value: 'P1' }, { label: 'N1', value: 'N1' }, { label: 'C', value: 'C' }]`
 */
export const getControlLabelOptions = (plate: Plate, measurement: string | null): PlateControlLabelOption[] => {
  const measurementStats = getMeasurementStats(plate, measurement)
  if (!measurementStats) return []

  return Object.keys(measurementStats).map((labelKey) => ({
    label: labelKey,
    value: labelKey,
  }))
}

/**
 * Picks legacy-compatible default positive/negative control labels.
 *
 * Priority:
 * - positive: `P` -> `P1` -> first available label
 * - negative: `N` -> `N1` -> second available label
 *
 * Returned data example:
 * - `{ positive: 'P1', negative: 'N1' }`
 */
export const getDefaultControlSelection = (
  plate: Plate,
  measurement: string | null,
): { positive: string | null; negative: string | null } => {
  const options = getControlLabelOptions(plate, measurement)
  const labels = options.map((option) => option.value)

  if (labels.length === 0) {
    return { positive: 'P', negative: 'N' }
  }

  const positive = labels.includes('P') ? 'P' : labels.includes('P1') ? 'P1' : (labels[0] ?? null)
  const negative = labels.includes('N') ? 'N' : labels.includes('N1') ? 'N1' : (labels[1] ?? null)

  return { positive, negative }
}

/**
 * Computes z-prime score from selected measurement and controls.
 *
 * Returned data example:
 * - `0.61`
 * - `null` when required stats are missing
 */
export const computeZPrime = (
  plate: Plate,
  measurement: string | null,
  positiveControl: string | null,
  negativeControl: string | null,
  timestampIndex: number,
): number | null => {
  const positiveStats = getPlateStatsForLabel(plate, measurement, positiveControl)
  const negativeStats = getPlateStatsForLabel(plate, measurement, negativeControl)
  if (!positiveStats || !negativeStats) return null

  const madPos = getStatsSeriesValue(positiveStats.mad, timestampIndex)
  const madNeg = getStatsSeriesValue(negativeStats.mad, timestampIndex)
  const medianPos = getStatsSeriesValue(positiveStats.median, timestampIndex)
  const medianNeg = getStatsSeriesValue(negativeStats.median, timestampIndex)

  if (madPos === null || madNeg === null || medianPos === null || medianNeg === null) return null
  const denominator = Math.abs(medianPos - medianNeg)
  if (denominator === 0) return null

  return 1 - (3 * (madPos + madNeg)) / denominator
}

/**
 * Computes SSMD score from selected measurement and controls.
 *
 * Returned data example:
 * - `44.18`
 * - `null` when required stats are missing
 */
export const computeSSMD = (
  plate: Plate,
  measurement: string | null,
  positiveControl: string | null,
  negativeControl: string | null,
  timestampIndex: number,
): number | null => {
  const positiveStats = getPlateStatsForLabel(plate, measurement, positiveControl)
  const negativeStats = getPlateStatsForLabel(plate, measurement, negativeControl)
  if (!positiveStats || !negativeStats) return null

  const madPos = getStatsSeriesValue(positiveStats.mad, timestampIndex)
  const madNeg = getStatsSeriesValue(negativeStats.mad, timestampIndex)
  const medianPos = getStatsSeriesValue(positiveStats.median, timestampIndex)
  const medianNeg = getStatsSeriesValue(negativeStats.median, timestampIndex)

  if (madPos === null || madNeg === null || medianPos === null || medianNeg === null) return null

  const denominator = 0.5 * (madPos + madNeg)
  if (denominator === 0) return null

  return Math.abs(medianPos - medianNeg) / denominator
}

/**
 * Reads overall min/max values for selected measurement at one timestamp.
 *
 * Returned data example:
 * - `{ min: 8, max: 1000 }`
 */
export const getOverallMinMaxForSelection = (
  plate: Plate,
  measurement: string | null,
  timestampIndex: number,
): { min: number; max: number } => {
  if (!measurement) {
    return { min: 0, max: 0 }
  }

  const overallStats = plate.details.overall_stats[measurement]
  if (!overallStats) {
    return { min: 0, max: 0 }
  }

  return {
    min: getStatsSeriesValue(overallStats.min, timestampIndex) ?? 0,
    max: getStatsSeriesValue(overallStats.max, timestampIndex) ?? 0,
  }
}
