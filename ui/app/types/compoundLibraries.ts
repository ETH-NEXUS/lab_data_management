import type { CompoundLibrary } from '~/types/lab'

/**
 * API constants and payloads for compound-library navigation features.
 *
 * Data examples:
 * - endpoint response: `{ results: [{ id: 5, name: 'Kinase Library', file_name: 'kinase.csv', plates: [] }] }`
 * - create plate payload: `{ barcode: 'LIB-001', library: 5 }`
 */

export const COMPOUND_LIBRARIES_ENDPOINT = 'compoundlibraries/'
export const COMPOUND_LIBRARY_ADD_PLATE_ENDPOINT = 'plates/'

export const COMPOUND_LIBRARIES_FETCH_ERROR_MESSAGE = 'Failed to load compound libraries.'
export const COMPOUND_LIBRARY_ADD_PLATE_ERROR_MESSAGE = 'Failed to create library plate.'

/**
 * Request payload for creating one plate inside one compound library.
 *
 * Data example:
 * - `{ barcode: 'LIB-001', library: 5 }`
 */
export type AddCompoundLibraryPlatePayload = {
  barcode: string
  library: number
}

/**
 * Helper function return value used after inserting a new plate into local state.
 *
 * Data example:
 * - `{ updatedLibrary: { id: 5, name: 'Kinase Library', file_name: 'kinase.csv', plates: [{ id: 41, barcode: 'LIB-001', ... }] }, createdPlateId: 41 }`
 */
export type AddCompoundLibraryPlateResult = {
  updatedLibrary: CompoundLibrary
  createdPlateId: number
}
