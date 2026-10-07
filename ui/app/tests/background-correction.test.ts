/**
 * Tests for finding the corrected measurements of a plate: the second heatmap of
 * the plate page shows them for the measurement selected in the main heatmap.
 */

import { describe, expect, it } from 'vitest'

import type { Plate } from '~/types/lab'
import { findBackgroundCorrections } from '~/utils/backgroundCorrection'

const STATS = { min: [0], max: [1], mean: [0], median: [0], std: [0], mad: [0] }

// Only the plate details are read
const plateWithLabels = (labels: string[], wellTypes: string[]): Plate => {
  const statsPerType = Object.fromEntries(wellTypes.map((type) => [type, STATS]))
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
      { label: 'Lum_bc_N1_median', referenceType: 'N1', method: 'median' },
      { label: 'Lum_bc_Nref_mean', referenceType: 'Nref', method: 'mean' },
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
})
