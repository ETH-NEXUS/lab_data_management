/**
 * API constants and payloads of the statistical analysis of an experiment.
 * Its output is read from the same long polling endpoint as a management command.
 *
 * Data examples:
 * - start payload: `{ experiment_id: 105, label: 'Lum_CTG', analysis_type: 'single', settings: {}, room_name: '3_1727600000000' }`
 * - results response: `{ results: ['20260929-125825_single_Lum_CTG.zip'] }`
 */

export const ANALYSIS_START_ENDPOINT = 'analysis/start/'
export const ANALYSIS_RESULTS_ENDPOINT = 'analysis/results/'
export const ANALYSIS_DOWNLOAD_ENDPOINT = 'analysis/download/'
export const ANALYSIS_CONDITIONS_ENDPOINT = 'analysis/conditions/'

export const ANALYSIS_START_ERROR_MESSAGE = 'Failed to start the analysis.'
export const ANALYSIS_RESULTS_ERROR_MESSAGE = 'Failed to load the analysis results.'
export const ANALYSIS_DOWNLOAD_ERROR_MESSAGE = 'Failed to download the analysis result.'
export const ANALYSIS_CONDITIONS_ERROR_MESSAGE = 'Failed to load the conditions of the measurement.'

export type AnalysisType = 'single' | 'selectivity'

/**
 * The conditions a selectivity analysis compares; empty for a single analysis.
 *
 * Data example:
 * - `{ condi_yes: 'irradiated', condi_no: 'not irradiated' }`
 */
export type AnalysisSettings = {
  condi_yes?: string
  condi_no?: string
}

/**
 * `positive_control` and `negative_control` are well types of the experiment,
 * e.g. 'P1' and 'N1'; the R report gets them as 'P' and 'N'.
 */
export type StartAnalysisPayload = {
  experiment_id: number
  label: string
  analysis_type: AnalysisType
  settings: AnalysisSettings
  positive_control: string
  negative_control: string
  room_name: string
}

export type AnalysisResultsResponse = {
  results: string[]
}

/**
 * Data example: `{ conditions: ['irradiated', 'not irradiated'] }`
 */
export type AnalysisConditionsResponse = {
  conditions: string[]
}
