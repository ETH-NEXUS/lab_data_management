/**
 * Tests for the colour scale of the heatmaps: the full range of the values, or
 * a robust range around the median that outliers do not stretch.
 */

import { describe, expect, it } from 'vitest'

import type { PlateStats } from '~/types/lab'
import { getHeatmapRange } from '~/utils/heatmapScale'

// The plate of the lab: values from 33 to 126340, median 734, MAD 400
const stats: PlateStats = { min: [33], max: [126340], mean: [4822], median: [734], std: [20000], mad: [400] }

describe('getHeatmapRange', () => {
  it('spans all values for the full range', () => {
    expect(getHeatmapRange({ Lum1: stats }, 'Lum1', 0, 'full')).toEqual({
      min: 33,
      max: 126340,
      lowerClipped: false,
      upperClipped: false,
    })
  })

  it('spans median ± 3 MAD for the robust range, within the values of the plate', () => {
    // 734 - 1200 is below the lowest value, so the lower end stays at 33
    expect(getHeatmapRange({ Lum1: stats }, 'Lum1', 0, 'robust')).toEqual({
      min: 33,
      max: 1934,
      lowerClipped: false,
      upperClipped: true,
    })
  })

  it('uses the full range when the wells do not spread', () => {
    const alike = { ...stats, mad: [0] }
    expect(getHeatmapRange({ Lum1: alike }, 'Lum1', 0, 'robust')).toEqual({
      min: 33,
      max: 126340,
      lowerClipped: false,
      upperClipped: false,
    })
  })

  it('has no range without statistics or without a measurement', () => {
    expect(getHeatmapRange({ Lum1: stats }, null, 0, 'full')).toEqual({
      min: 0,
      max: 0,
      lowerClipped: false,
      upperClipped: false,
    })
    expect(getHeatmapRange(undefined, 'Lum1', 0, 'robust')).toEqual({
      min: 0,
      max: 0,
      lowerClipped: false,
      upperClipped: false,
    })
  })
})
