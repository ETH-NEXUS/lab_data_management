/**
 * API constants and payloads of the background correction of one plate.
 * The corrected values are saved as a new measurement of the plate, named
 * `<label>_bc_<reference type>_<method>`, e.g. `Lum_CTG_bc_N1_median`.
 *
 * Data examples:
 * - request: `{ label: 'Lum_CTG', reference_type: 'N1', method: 'median' }`
 * - response: `{ label: 'Lum_CTG_bc_N1_median' }`
 */

export const BACKGROUND_CORRECTION_ENDPOINT = 'background_correction/plates/'

export const BACKGROUND_CORRECTION_ERROR_MESSAGE = 'Failed to correct the background of the plate.'

export const BACKGROUND_CORRECTION_METHODS = ['median', 'mean'] as const

export type BackgroundCorrectionMethod = (typeof BACKGROUND_CORRECTION_METHODS)[number]

export type BackgroundCorrectionSettings = {
  label: string
  reference_type: string
  method: BackgroundCorrectionMethod
}

export type BackgroundCorrectionResponse = {
  label: string
}

/**
 * One corrected measurement found on a plate.
 *
 * Data example:
 * - `{ label: 'Lum_CTG_bc_N1_median', referenceType: 'N1', method: 'median' }`
 */
export type BackgroundCorrection = {
  label: string
  referenceType: string
  method: BackgroundCorrectionMethod
}

/**
 * The log10 of a measurement of one plate, saved as `<label>_log10`.
 * Wells with a value of 0 or below are left empty and counted.
 *
 * Data examples:
 * - request: `{ label: 'Lum1' }`
 * - response: `{ label: 'Lum1_log10', skipped: 2 }`
 */
export const LOG10_ENDPOINT = 'log10/'
export const LOG10_ERROR_MESSAGE = 'Failed to calculate the log10 of the measurement.'
export const LOG10_SUFFIX = '_log10'

export type Log10Response = {
  label: string
  skipped: number
}

// The calculations the plate page offers
export const PLATE_CALCULATIONS = ['background_correction', 'log10'] as const

export type PlateCalculation = (typeof PLATE_CALCULATIONS)[number]

/**
 * What the calculation window reports after a calculation.
 *
 * Data example:
 * - `{ calculation: 'log10', label: 'Lum1', newLabel: 'Lum1_log10' }`
 */
export type PlateCalculationResult = {
  calculation: PlateCalculation
  label: string
  newLabel: string
}
