import type { Plate, PlateLabelValue } from '~/types/lab'

/**
 * Shared plate-page UI constants and lightweight page-level types.
 *
 * Data examples:
 * - `PLATE_PAGE_DEFAULT_LEFT_PERCENT = 58`
 * - `PLATE_PAGE_MIN_LEFT_PERCENT = 28`
 * - `PLATE_PAGE_MAX_LEFT_PERCENT = 78`
 */
export const PLATE_PAGE_DEFAULT_LEFT_PERCENT = 58
export const PLATE_PAGE_MIN_LEFT_PERCENT = 28
export const PLATE_PAGE_MAX_LEFT_PERCENT = 78

/**
 * Reusable props for plate pane heading blocks.
 *
 * Accepted data example:
 * - `{ title: 'demo_1', description: 'Dynamic plate visualization will be migrated here next.' }`
 */
export type PlatePaneHeaderProps = {
  title: string
  description: string
  titleClass?: string
  descriptionClass?: string
}

/**
 * Plate barcode accepted by the migrated plate page/store.
 *
 * Accepted data examples:
 * - `'demo_1'`
 * - `'250401control'`
 */
export type PlateBarcode = string

/**
 * API endpoint constants used by the plate store.
 */
export const PLATES_ENDPOINT = 'plates/'
export const PLATE_BARCODES_ENDPOINT = 'plates/barcodes/'

/**
 * Store-level fallback messages for plate page fetch flows.
 */
export const PLATE_FETCH_ERROR_MESSAGE = 'Failed to load plate.'
export const PLATE_TARGET_BARCODES_FETCH_ERROR_MESSAGE = 'Failed to load target plate barcode options.'
export const PLATE_TEMPLATE_BARCODES_FETCH_ERROR_MESSAGE = 'Failed to load template plate barcode options.'

/**
 * Builds a not-found message for one barcode.
 *
 * Returned message example:
 * - `'No plate found with barcode demo_1.'`
 */
export const getPlateNotFoundMessage = (barcode: string): string => `No plate found with barcode ${barcode}.`

/**
 * Builds a duplicated-result message when a barcode query matches more than one row.
 *
 * Returned message example:
 * - `'Multiple plates found with barcode demo_1.'`
 */
export const getMultiplePlatesFoundMessage = (barcode: string): string =>
  `Multiple plates found with barcode ${barcode}.`

/**
 * Aggregated plate-page fetch payload returned by `initializePlatePage`.
 *
 * Returned data example:
 * - `{ plate: { id: 21, barcode: 'demo_1' }, targetPlateBarcodeOptions: [{ label: 'A001', value: 11 }], templatePlateBarcodeOptions: [] }`
 */
export type PlatePageData = {
  plate: Plate | null
  targetPlateBarcodeOptions: PlateLabelValue[]
  templatePlateBarcodeOptions: PlateLabelValue[]
}
