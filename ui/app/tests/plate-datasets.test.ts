/**
 * Tests for finding how the measurements of a plate were calculated (background
 * correction, log10, normalization), and for the compact statistics below the heatmap.
 */

import { describe, expect, it } from 'vitest'

import type { Plate } from '~/types/lab'
import {
  countWellsWithoutLog10,
  findBackgroundCorrections,
  formatSummaryNumber,
  getCorrectionDataset,
  getLog10SourceLabel,
  getNormalizedDataset,
  getWellTypeLog10Median,
  getWellTypeStatistic,
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
    // One read of each measurement
    measurement_timestamps: { Lum1: ['2026-09-30T15:48:08'], Lum1_log10: ['2026-09-30T15:48:08'] },
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

  it('counts the wells of one type left empty by the log10', () => {
    expect(countWellsWithoutLog10(plate, 'Lum1', 'Lum1_log10', 'N')).toBe(1)
    expect(countWellsWithoutLog10(plate, 'Lum1', 'Lum1_log10', 'R')).toBe(0)
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

  it('does not take a statistic of a well type that misses a read of the plate', () => {
    // Two reads of the plate, but the R wells only have a value in one of them
    const twoReads = {
      ...plate,
      details: {
        ...plate.details,
        measurement_timestamps: { Lum1: ['2026-09-30T15:48:08', '2026-09-30T17:00:00'] },
        stats: { Lum1: { N: { ...stats(2958), median: [2958, 3100] }, R: stats(734) } },
      },
    } as unknown as Plate

    expect(getWellTypeStatistic(twoReads, 'Lum1', 'N', 'median', 1)).toBe(3100)
    expect(getWellTypeStatistic(twoReads, 'Lum1', 'R', 'median', 0)).toBeNull()
  })

  it('finds how a %Inhibition or %Activity was calculated', () => {
    expect(getNormalizedDataset(plate, 'Lum1_activity_N_P')).toEqual({
      kind: 'activity',
      label: 'Lum1_activity_N_P',
      source: 'Lum1',
      negativeType: 'N',
      positiveType: 'P',
      inhibitionLabel: 'Lum1_inhibition_N_P',
      activityLabel: 'Lum1_activity_N_P',
    })
    expect(getNormalizedDataset(plate, 'Lum1_inhibition_N_P')?.kind).toBe('inhibition')
    expect(getNormalizedDataset(plate, 'Lum1')).toBeNull()
  })

  it('takes the median of the log10 values of a well type, without wells of 0 or below', () => {
    // The N wells have the raw values 0 and 1000: the well with 0 has no log10
    expect(getWellTypeLog10Median(plate, 'Lum1', 'N', 0)).toBe(3)
  })
})

// Only the plate details are read
const plateWithLabels = (labels: string[], wellTypes: string[]): Plate => {
  const statsPerType = Object.fromEntries(wellTypes.map((type) => [type, stats(0)]))
  return {
    details: {
      measurement_labels: labels,
      measurement_timestamps: {},
      stats: { Lum: statsPerType },
      overall_stats: {},
    },
  } as unknown as Plate
}

describe('findBackgroundCorrections', () => {
  it('finds every reference and method of the measurement', () => {
    const plate = plateWithLabels(['Lum', 'Lum_bc_N1_median', 'Lum_bc_Nref_mean'], ['C', 'N1', 'Nref', 'P'])

    expect(findBackgroundCorrections(plate, 'Lum')).toEqual([
      { label: 'Lum_bc_N1_median', source: 'Lum', referenceType: 'N1', method: 'median' },
      { label: 'Lum_bc_Nref_mean', source: 'Lum', referenceType: 'Nref', method: 'mean' },
    ])
  })

  it('does not take the corrections of another measurement', () => {
    const plate = plateWithLabels(['Lum', 'Lum2_bc_N1_median', 'Fluo_bc_N1_median'], ['C', 'N1'])

    expect(findBackgroundCorrections(plate, 'Lum')).toEqual([])
  })

  it('finds nothing without a selected measurement', () => {
    const plate = plateWithLabels(['Lum', 'Lum_bc_N1_median'], ['C', 'N1'])

    expect(findBackgroundCorrections(plate, null)).toEqual([])
  })

  it('finds how a correction shown in the main heatmap was calculated', () => {
    const plate = plateWithLabels(['Lum', 'Lum_bc_N1_median'], ['C', 'N1'])

    expect(getCorrectionDataset(plate, 'Lum_bc_N1_median')).toEqual({
      label: 'Lum_bc_N1_median',
      source: 'Lum',
      referenceType: 'N1',
      method: 'median',
    })
    expect(getCorrectionDataset(plate, 'Lum')).toBeNull()
  })
})
