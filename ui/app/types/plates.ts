import type { Plate, PlateLabelValue, PlateStats, WellInfo } from '~/types/lab'

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
 * Allowed modes for what is shown inside one plate well.
 *
 * Data examples:
 * - `'hr_position'`
 * - `'type'`
 * - `'position'`
 */
export type PlateWellContentMode = 'hr_position' | 'type' | 'position'

/**
 * Heatmap palette value (from/to color).
 *
 * Data example:
 * - `{ from: '#FF0000', to: '#00FF00' }`
 */
export type PlatePaletteValue = {
  from: string
  to: string
}

/**
 * Selectable heatmap palette option.
 *
 * Data example:
 * - `{ label: 'GreenRed', value: { from: '#FF0000', to: '#00FF00' } }`
 */
export type PlatePaletteOption = {
  label: string
  value: PlatePaletteValue
}

/**
 * UI option for positive/negative control selectors.
 *
 * Data example:
 * - `{ label: 'P1', value: 'P1' }`
 */
export type PlateControlLabelOption = {
  label: string
  value: string
}

/**
 * UI option for timestamp selector.
 *
 * Data example:
 * - `{ label: '2026-01-01T00:00:00Z', value: 0 }`
 */
export type PlateTimestampOption = {
  label: string
  value: number
}

/**
 * Minimal state handled by the dedicated plate-view UI store.
 *
 * Data example:
 * - `{ splitter: 58, selectedWellInfo: undefined, wellContent: 'hr_position', showHeatmap: false, smallerMapView: false, heatmapPalette: { label: 'GreenRed', value: { from: '#FF0000', to: '#00FF00' } }, squareCompoundType: false }`
 */
export type PlateViewState = {
  splitter: number
  selectedWellInfo: WellInfo | undefined
  measurementOptions: string[] | null
  selectedMeasurement: string | null
  selectedTimestampIdx: number
  selectedPosControl: string | null
  selectedNegControl: string | null
  wellContent: PlateWellContentMode
  showHeatmap: boolean
  smallerMapView: boolean
  heatmapPalette: PlatePaletteOption
  perPlateView: boolean
  squareCompoundType: boolean
  plotView: boolean
  showStructure: boolean
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

/**
 * Mapping of control label -> stats arrays for one selected measurement.
 *
 * Data example:
 * - `{ P1: { min: [1], max: [5], mean: [3], median: [3], std: [0.2], mad: [0.1] }, N1: { min: [0], max: [2], mean: [1], median: [1], std: [0.2], mad: [0.1] } }`
 */
export type PlateMeasurementStatsMap = Record<string, PlateStats>
