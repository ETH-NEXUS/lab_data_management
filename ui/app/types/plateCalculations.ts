/**
 * API constants and data of the calculations with the measurements of one plate
 * (plate page, "Calculations"). Each result is saved as a new measurement of the
 * plate, named after how it was calculated:
 * - background correction: `<label>_bc_<reference type>_<method>`, e.g. `Lum1_bc_R_median`
 * - log10: `<label>_log10`, e.g. `Lum1_log10`
 * - %Activity: `<label>_activity_<negative>_<positive>`, e.g. `Lum1_activity_N_P`
 */

export const PLATE_CALCULATIONS_ENDPOINT = 'plate_calculations/plates/'
export const PLATE_CALCULATION_ERROR_MESSAGE = 'The calculation of the plate failed.'

// The calculations the plate page offers
export const PLATE_CALCULATIONS = ['background_correction', 'log10', 'percent_activity'] as const

export type PlateCalculation = (typeof PLATE_CALCULATIONS)[number]

// The last part of the endpoint of each calculation, e.g. `plate_calculations/plates/42/log10/`
export const PLATE_CALCULATION_PATHS: Record<PlateCalculation, string> = {
  background_correction: 'background_correction/',
  log10: 'log10/',
  percent_activity: 'activity/',
}

// The parts of the names of the new measurements
export const CORRECTION_INFIX = '_bc_'
export const LOG10_SUFFIX = '_log10'
export const ACTIVITY_INFIX = '_activity_'

export const BACKGROUND_CORRECTION_METHODS = ['median', 'mean'] as const

export type BackgroundCorrectionMethod = (typeof BACKGROUND_CORRECTION_METHODS)[number]

/**
 * What is sent to the server for each calculation.
 *
 * Data examples:
 * - background correction: `{ label: 'Lum1', reference_type: 'R', method: 'median' }`
 * - log10: `{ label: 'Lum1' }`
 * - %Activity: `{ label: 'Lum1', negative_type: 'N', positive_type: 'P' }`
 */
export type PlateCalculationSettings =
  | { label: string; reference_type: string; method: BackgroundCorrectionMethod }
  | { label: string }
  | { label: string; negative_type: string; positive_type: string }

/**
 * The answer of the server: the name of the new measurement, and for log10 how
 * many wells were left empty because a value is 0 or below.
 *
 * Data examples:
 * - `{ label: 'Lum1_bc_R_median' }`
 * - `{ label: 'Lum1_log10', skipped: 2 }`
 */
export type PlateCalculationResponse = {
  label: string
  skipped?: number
}

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

/**
 * A background correction found on a plate.
 *
 * Data example:
 * - `{ label: 'Lum1_bc_R_median', source: 'Lum1', referenceType: 'R', method: 'median' }`
 */
export type BackgroundCorrection = {
  label: string
  source: string
  referenceType: string
  method: BackgroundCorrectionMethod
}

/**
 * A %Activity found on a plate.
 *
 * Data example:
 * - `{ label: 'Lum1_activity_N_P', source: 'Lum1', negativeType: 'N', positiveType: 'P' }`
 */
export type ActivityDataset = {
  label: string
  source: string
  negativeType: string
  positiveType: string
}
