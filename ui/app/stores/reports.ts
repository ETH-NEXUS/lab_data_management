import { defineStore } from 'pinia'
import { ref } from 'vue'
import {
  REPORT_DOWNLOAD_ENDPOINT,
  REPORT_DOWNLOAD_ERROR_MESSAGE,
  REPORT_GENERATE_ENDPOINT,
  REPORT_GENERATE_ERROR_MESSAGE,
  REPORT_INPUT_NOTEBOOKS_DIR,
  REPORT_INPUT_NOTEBOOKS_ERROR_MESSAGE,
  REPORT_INPUT_NOTEBOOKS_FILE_FORMAT,
  REPORT_LIST_FILES_ENDPOINT,
  REPORT_OUTPUT_NOTEBOOKS_ERROR_MESSAGE,
  type GenerateReportPayload,
  type ListInputNotebooksPayload,
  type ListOutputReportsPayload,
  type ReportFilesResponse,
} from '~/types/reports'
import { requestApiData } from '~/utils/apiRequests'

export const useReportStore = defineStore('reportStore', () => {
  const inputNotebookOptions = ref<string[]>([])
  const outputNotebooks = ref<string[]>([])
  const isLoadingInputNotebookOptions = ref(false)
  const isLoadingOutputNotebooks = ref(false)
  const isGeneratingReport = ref(false)
  const isDownloadingReport = ref(false)

  /**
   * Extracts a filename from a full notebook/report path.
   *
   * Accepted input examples:
   * - `'/notebooks/input/template.ipynb'`
   * - `'/notebooks/output/exp_1.pdf'`
   *
   * Returned data examples:
   * - `'template.ipynb'`
   * - `'exp_1.pdf'`
   */
  const getFileNameFromPath = (path: string): string => {
    const pathParts = path.split('/')
    if (pathParts.length === 0) {
      return path
    }

    const lastPart = pathParts[pathParts.length - 1]
    if (!lastPart) {
      return path
    }

    return lastPart
  }

  /**
   * Loads available notebook templates from `/notebooks/input`.
   *
   * Returned data example:
   * - `['/notebooks/input/default_report.ipynb']`
   */
  const fetchInputNotebookOptions = async (): Promise<string[]> => {
    isLoadingInputNotebookOptions.value = true
    try {
      const payload: ListInputNotebooksPayload = {
        notebooks_dir: REPORT_INPUT_NOTEBOOKS_DIR,
        file_format: REPORT_INPUT_NOTEBOOKS_FILE_FORMAT,
      }

      const response = await requestApiData<ReportFilesResponse>(
        REPORT_LIST_FILES_ENDPOINT,
        {
          method: 'POST',
          body: payload,
        },
        REPORT_INPUT_NOTEBOOKS_ERROR_MESSAGE,
      )

      inputNotebookOptions.value = response.notebooks ?? []
      return inputNotebookOptions.value
    } finally {
      isLoadingInputNotebookOptions.value = false
    }
  }

  /**
   * Loads generated report files for one experiment.
   *
   * Accepted input example:
   * - `'Dose response'`
   *
   * Returned data example:
   * - `['/notebooks/output/Dose response_report.pdf']`
   */
  const fetchOutputNotebooks = async (experimentName: string): Promise<string[]> => {
    const trimmedExperimentName = experimentName.trim()
    if (trimmedExperimentName === '') {
      outputNotebooks.value = []
      return outputNotebooks.value
    }

    isLoadingOutputNotebooks.value = true
    try {
      const payload: ListOutputReportsPayload = {
        experiment: trimmedExperimentName,
      }

      const response = await requestApiData<ReportFilesResponse>(
        REPORT_LIST_FILES_ENDPOINT,
        {
          method: 'POST',
          body: payload,
        },
        REPORT_OUTPUT_NOTEBOOKS_ERROR_MESSAGE,
      )

      outputNotebooks.value = response.notebooks ?? []
      return outputNotebooks.value
    } finally {
      isLoadingOutputNotebooks.value = false
    }
  }

  /**
   * Triggers PDF report generation and refreshes the output report list.
   *
   * Accepted payload example:
   * - `{ experiment: 'Dose response', label: 'OD600', notebook_path: '/notebooks/input/default.ipynb', selected_pos: 'P', selected_neg: 'N' }`
   */
  const generateReport = async (payload: GenerateReportPayload): Promise<void> => {
    isGeneratingReport.value = true
    try {
      await requestApiData<unknown>(
        REPORT_GENERATE_ENDPOINT,
        {
          method: 'POST',
          body: payload,
        },
        REPORT_GENERATE_ERROR_MESSAGE,
      )

      await fetchOutputNotebooks(payload.experiment)
    } finally {
      isGeneratingReport.value = false
    }
  }

  /**
   * Downloads one generated PDF report.
   *
   * Accepted input example:
   * - `'/notebooks/output/Dose response_report.pdf'`
   */
  const downloadReport = async (path: string): Promise<void> => {
    isDownloadingReport.value = true
    try {
      if (typeof window === 'undefined' || typeof document === 'undefined') {
        return
      }

      const reportBlob = await requestApiData<Blob>(
        REPORT_DOWNLOAD_ENDPOINT,
        {
          method: 'POST',
          body: {
            path,
          },
          responseType: 'blob',
        },
        REPORT_DOWNLOAD_ERROR_MESSAGE,
      )

      const downloadFileName = getFileNameFromPath(path)
      const reportBlobUrl = window.URL.createObjectURL(reportBlob)
      const link = document.createElement('a')
      link.href = reportBlobUrl
      link.download = downloadFileName
      document.body.appendChild(link)
      link.click()
      document.body.removeChild(link)
      window.URL.revokeObjectURL(reportBlobUrl)
    } finally {
      isDownloadingReport.value = false
    }
  }

  return {
    inputNotebookOptions,
    outputNotebooks,
    isLoadingInputNotebookOptions,
    isLoadingOutputNotebooks,
    isGeneratingReport,
    isDownloadingReport,
    fetchInputNotebookOptions,
    fetchOutputNotebooks,
    generateReport,
    downloadReport,
    getFileNameFromPath,
  }
})
