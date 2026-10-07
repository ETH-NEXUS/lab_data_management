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
