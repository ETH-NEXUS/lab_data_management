import type { Plate } from '~/types/lab'
import { BACKGROUND_CORRECTION_METHODS, type BackgroundCorrection } from '~/types/backgroundCorrection'

/**
 * Name of the corrected measurement, the same as the server gives it.
 *
 * Example: `('Lum_CTG', 'N1', 'median')` -> `'Lum_CTG_bc_N1_median'`
 */
export const getCorrectedLabel = (label: string, referenceType: string, method: string): string => {
  return `${label}_bc_${referenceType}_${method}`
}

/**
 * Well types that have a value of the measurement on the plate, e.g. `['C', 'N1', 'P']`.
 */
export const getWellTypesOfMeasurement = (plate: Plate, label: string | null): string[] => {
  if (!label) return []
  return Object.keys(plate.details.stats[label] ?? {}).sort()
}

/**
 * Corrected measurements of one measurement of the plate.
 * Each well type of the measurement could have been the reference, so every
 * well type and method is looked up by its corrected label.
 *
 * Returned data example (plate with `Lum_CTG_bc_N1_median` and `Lum_CTG_bc_N1_mean`):
 * - `[{ label: 'Lum_CTG_bc_N1_median', referenceType: 'N1', method: 'median' },
 *     { label: 'Lum_CTG_bc_N1_mean', referenceType: 'N1', method: 'mean' }]`
 */
export const findBackgroundCorrections = (plate: Plate, label: string | null): BackgroundCorrection[] => {
  if (!label) return []

  const plateLabels = plate.details.measurement_labels ?? []
  const corrections: BackgroundCorrection[] = []

  for (const referenceType of getWellTypesOfMeasurement(plate, label)) {
    for (const method of BACKGROUND_CORRECTION_METHODS) {
      const correctedLabel = getCorrectedLabel(label, referenceType, method)
      if (plateLabels.includes(correctedLabel)) {
        corrections.push({ label: correctedLabel, referenceType, method })
      }
    }
  }
  return corrections
}
