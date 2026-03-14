import type { BarcodeSpecification, ExperimentDetails, Plate } from '~/types/lab'

/**
 * Minimal experiment model used by list/detail queries and store updates.
 *
 * Data example:
 * - `{ id: 21, name: 'Dose response', project: 5, available_measurement_labels: ['OD600'], details: { d: 384, project_id: 5, measurement_labels: ['OD600'], measurement_timestamps: { OD600: ['2026-01-01T00:00:00Z'] }, stats: {}, overall_stats: {} }, plates: [] }`
 */
export type Experiment = {
  id: number
  name: string
  project: number
  description?: string | null
  created_at?: Date | string | null
  plates: Plate[]
  available_measurement_labels: string[]
  details: ExperimentDetails
  barcode_specifications?: BarcodeSpecification[]
  [key: string]: unknown
}

export const EXPERIMENTS_QUERY_KEY = ['experiments'] as const
export const EXPERIMENTS_ENDPOINT = 'experiments/'
export const EXPERIMENTS_ERROR_MESSAGE = 'Failed to load experiments.'
export const EXPERIMENT_ERROR_MESSAGE = 'Failed to load experiment.'
export const EXPERIMENT_CREATE_ERROR_MESSAGE = 'Failed to create experiment.'
export const EXPERIMENT_UPDATE_ERROR_MESSAGE = 'Failed to update experiment.'
export const EXPERIMENT_BARCODE_CREATE_ERROR_MESSAGE = 'Failed to create barcode specification.'
export const EXPERIMENT_BARCODE_DELETE_ERROR_MESSAGE = 'Failed to delete barcode specification.'

/**
 * Query key helper for one experiment.
 *
 * Data example:
 * - `getExperimentQueryKey(21)` returns `['experiment', '21']`
 */
export const getExperimentQueryKey = (experimentId: number | string) => ['experiment', String(experimentId)] as const

/**
 * PATCH payload for experiment updates.
 *
 * Accepted payload examples:
 * - `{ name: 'Updated experiment name' }`
 * - `{ description: 'Updated experiment description' }`
 */
export type ExperimentPayload = {
  name?: string
  description?: string | null
  [key: string]: unknown
}

/**
 * Request payload for creating an experiment.
 *
 * Accepted payload examples:
 * - `{ name: 'Dose response', project: 5 }`
 * - `{ name: 'QC run', project: 8 }`
 */
export type CreateExperimentPayload = {
  name: string
  project: number
}

/**
 * Allowed plate sides used for barcode generation.
 *
 * Data examples:
 * - `'North'`
 * - `'West'`
 */
export type BarcodeSide = 'North' | 'South' | 'East' | 'West'

/**
 * Request payload for creating a barcode specification.
 *
 * Accepted payload examples:
 * - `{ experiment_id: 21, prefix: 'A001', number_of_plates: 12, sides: ['North', 'South'] }`
 * - `{ experiment_id: 22, prefix: 'B', number_of_plates: 4, sides: ['East'] }`
 */
export type CreateBarcodeSpecificationPayload = {
  experiment_id: number
  prefix: string
  number_of_plates: number
  sides: BarcodeSide[]
}
