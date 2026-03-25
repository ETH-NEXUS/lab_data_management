/**
 * Shared report API constants and payload types.
 *
 * Data structure examples:
 * - `{ notebooks: ['/notebooks/input/template.ipynb'] }`
 * - `{ experiment: 'Dose response' }`
 * - `{ path: '/notebooks/output/report_1.pdf' }`
 */

export const REPORT_LIST_FILES_ENDPOINT = 'list_files/'
export const REPORT_GENERATE_ENDPOINT = 'generate_pdf_report/'
export const REPORT_DOWNLOAD_ENDPOINT = 'download_pdf_report/'

export const REPORT_INPUT_NOTEBOOKS_DIR = '/notebooks/input'
export const REPORT_INPUT_NOTEBOOKS_FILE_FORMAT = '.ipynb'

export const REPORT_INPUT_NOTEBOOKS_ERROR_MESSAGE = 'Failed to load notebook templates.'
export const REPORT_OUTPUT_NOTEBOOKS_ERROR_MESSAGE = 'Failed to load generated reports.'
export const REPORT_GENERATE_ERROR_MESSAGE = 'Failed to generate report.'
export const REPORT_DOWNLOAD_ERROR_MESSAGE = 'Failed to download report.'

/**
 * Response body returned by `POST /api/list_files/`.
 *
 * Data example:
 * - `{ notebooks: ['/notebooks/input/template.ipynb', '/notebooks/output/report.pdf'] }`
 */
export type ReportFilesResponse = {
  notebooks: string[]
}

/**
 * Request body for loading input notebook templates.
 *
 * Data example:
 * - `{ notebooks_dir: '/notebooks/input', file_format: '.ipynb' }`
 */
export type ListInputNotebooksPayload = {
  notebooks_dir: string
  file_format: string
}

/**
 * Request body for loading generated report files for one experiment.
 *
 * Data example:
 * - `{ experiment: 'Dose response' }`
 */
export type ListOutputReportsPayload = {
  experiment: string
}

/**
 * Request body for generating one report file.
 *
 * Data example:
 * - `{ experiment: 'Dose response', label: 'OD600', notebook_path: '/notebooks/input/report.ipynb', selected_pos: 'P', selected_neg: 'N' }`
 */
export type GenerateReportPayload = {
  experiment: string
  label: string
  notebook_path: string
  selected_pos: string
  selected_neg: string
}

/**
 * Request body for downloading one generated report.
 *
 * Data example:
 * - `{ path: '/notebooks/output/dose_response.pdf' }`
 */
export type DownloadReportPayload = {
  path: string
}
