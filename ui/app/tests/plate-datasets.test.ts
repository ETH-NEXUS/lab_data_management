/**
 * Tests for the compact statistics and the log10 report below the heatmap of a plate.
 */

import { describe, expect, it } from 'vitest'

import type { Plate } from '~/types/lab'
import {
  countWellsWithoutLog10,
  formatSummaryNumber,
  getActivityDataset,
  getLog10SourceLabel,
  getWellTypeMedian,
  summarizeDataset,
} from '~/utils/plateDatasets'

const stats = (median: number) => ({ min: [0], max: [0], mean: [0], median: [median], std: [0], mad: [0] })

// Three wells; the second one has a raw value of 0, so no log10
const plate = {
  wells: [
    { type: 'R', measurements: { Lum1: [100], Lum1_log10: [2] } },
    { type: 'N', measurements: { Lum1: [0] } },
    { type: 'N', measurements: { Lum1: [1000], Lum1_log10: [3] } },
  ],
  details: {
    measurement_labels: ['Lum1', 'Lum1_log10', 'Lum1_activity_N_P'],
    measurement_timestamps: {},
    stats: { Lum1: { N: stats(2958), P: stats(673) }, Lum1_log10: { R: stats(2), N: stats(3) } },
    overall_stats: { Lum1_log10: { min: [2], max: [3], mean: [2.5], median: [2.5], std: [0.5], mad: [0.5] } },
  },
} as unknown as Plate

describe('plate datasets', () => {
  it('finds the measurement a log10 was calculated from', () => {
    expect(getLog10SourceLabel(plate, 'Lum1_log10')).toBe('Lum1')
    expect(getLog10SourceLabel(plate, 'Lum1')).toBeNull()
    expect(getLog10SourceLabel(plate, 'Fluo_log10')).toBeNull()
  })

  it('counts the wells left empty by the log10', () => {
    expect(countWellsWithoutLog10(plate, 'Lum1', 'Lum1_log10')).toBe(1)
  })

  it('summarizes a measurement at one time point', () => {
    expect(summarizeDataset(plate, 'Lum1_log10', 0)).toEqual({
      wells: 2,
      median: 2.5,
      min: 2,
      max: 3,
      byWellType: [
        { wellType: 'N', median: 3 },
        { wellType: 'R', median: 2 },
      ],
    })
  })

  it('shows at most 3 decimals', () => {
    expect(formatSummaryNumber(2.86634)).toBe('2.866')
    expect(formatSummaryNumber(734)).toBe('734')
    expect(formatSummaryNumber(null)).toBe('–')
  })

  it('finds how a %Activity was calculated', () => {
    expect(getActivityDataset(plate, 'Lum1_activity_N_P')).toEqual({
      label: 'Lum1_activity_N_P',
      source: 'Lum1',
      negativeType: 'N',
      positiveType: 'P',
    })
    expect(getActivityDataset(plate, 'Lum1')).toBeNull()
    expect(getWellTypeMedian(plate, 'Lum1', 'N', 0)).toBe(2958)
  })
})
