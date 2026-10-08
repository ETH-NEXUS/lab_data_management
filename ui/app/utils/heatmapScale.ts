import type { PlateStats } from '~/types/lab'
import type { HeatmapScale } from '~/types/plates'
import { getStatsSeriesValue } from '~/utils/plateStats'

// How many MADs around the median the robust scale covers. The MAD from the
// server is scaled like a standard deviation, so ± 3 MAD holds the usual wells.
const ROBUST_MADS = 3

/**
 * The values at the two ends of the colour scale. A clipped end is inside the
 * values of the plate: wells beyond it get the end colour.
 *
 * Data example (robust, an outlier of 126340 above):
 * - `{ min: 1.2, max: 5.0, lowerClipped: false, upperClipped: true }`
 */
export type HeatmapRange = {
  min: number
  max: number
  lowerClipped: boolean
  upperClipped: boolean
}

/**
 * The range of the colour scale of one measurement at one time point, from the
 * statistics the server keeps of a plate or of all plates of an experiment
 * (`overall_stats`, by measurement).
 *
 * Example: min 33, max 126340, median 734, MAD 400, `robust` -> `{ min: 33, max: 1934, lowerClipped: false, upperClipped: true }`
 */
export const getHeatmapRange = (
  statsByLabel: Record<string, PlateStats> | undefined,
  label: string | null,
  timestampIndex: number,
  scale: HeatmapScale,
): HeatmapRange => {
  const stats = label ? statsByLabel?.[label] : undefined
  const min = getStatsSeriesValue(stats?.min, timestampIndex) ?? 0
  const max = getStatsSeriesValue(stats?.max, timestampIndex) ?? 0
  const fullRange = { min, max, lowerClipped: false, upperClipped: false }
  if (scale === 'full') return fullRange

  const median = getStatsSeriesValue(stats?.median, timestampIndex)
  const mad = getStatsSeriesValue(stats?.mad, timestampIndex)
  // Without a spread (e.g. all wells alike) the robust range would be empty
  if (median === null || mad === null || mad === 0) return fullRange

  // The scale never reaches beyond the values of the plate
  const lower = Math.max(min, median - ROBUST_MADS * mad)
  const upper = Math.min(max, median + ROBUST_MADS * mad)
  return { min: lower, max: upper, lowerClipped: lower > min, upperClipped: upper < max }
}
