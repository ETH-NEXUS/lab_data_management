import { defineStore } from 'pinia'
import { ref } from 'vue'
import {
  ANALYSIS_DOWNLOAD_ENDPOINT,
  ANALYSIS_DOWNLOAD_ERROR_MESSAGE,
  ANALYSIS_RESULTS_ENDPOINT,
  ANALYSIS_RESULTS_ERROR_MESSAGE,
  ANALYSIS_START_ENDPOINT,
  ANALYSIS_START_ERROR_MESSAGE,
  type AnalysisResultsResponse,
  type StartAnalysisPayload,
} from '~/types/analysis'
import {
  MANAGEMENT_LONG_POLLING_ENDPOINT,
  MANAGEMENT_LONG_POLLING_ERROR_MESSAGE,
  type CommandMessage,
  type CommandStatus,
  type LongPollingResponse,
} from '~/types/management'
import { requestApiData, requestApiVoid } from '~/utils/apiRequests'
import { getErrorMessage } from '~/utils/errors'

const POLL_INTERVAL_MS = 2000
// A failed request for the output is repeated this many times before giving up
const MAX_FAILED_OUTPUT_REQUESTS = 5

const wait = (milliseconds: number) => new Promise((resolve) => setTimeout(resolve, milliseconds))

/**
 * The statistical analysis of an experiment: it runs in the celery container, the
 * page shows its messages while it runs and lists the result zips afterwards.
 */
export const useAnalysisStore = defineStore('analysisStore', () => {
  // e.g. [{ level: 'info', text: 'Step 1 of 3: collecting the data ...' }]
  const messages = ref<CommandMessage[]>([])
  const status = ref<CommandStatus | null>(null)
  const isRunning = ref(false)
  // The experiment of the last run: its messages are only shown on that experiment
  const experimentId = ref<number | null>(null)
  // e.g. ['20260929-125825_single_Lum_CTG.zip'], newest first
  const results = ref<string[]>([])

  /**
   * Starts one analysis and shows its messages until it has ended.
   *
   * Accepted data example:
   * - `{ experiment_id: 105, label: 'Lum_CTG', analysis_type: 'single', settings: {} }`
   */
  const startAnalysis = async (payload: Omit<StartAnalysisPayload, 'room_name'>): Promise<void> => {
    // A new room for every run, so the output of an earlier run is never shown
    const roomName = `${payload.experiment_id}_${Date.now()}`
    experimentId.value = payload.experiment_id
    messages.value = []
    status.value = 'running'
    isRunning.value = true

    try {
      await requestApiVoid(
        ANALYSIS_START_ENDPOINT,
        { method: 'POST', body: { ...payload, room_name: roomName } },
        ANALYSIS_START_ERROR_MESSAGE,
      )
    } catch (err: unknown) {
      messages.value = [{ level: 'error', text: `${ANALYSIS_START_ERROR_MESSAGE} ${getErrorMessage(err)}` }]
      status.value = 'failed'
      isRunning.value = false
      return
    }

    await readOutput(roomName)
  }

  /**
   * Adds the new messages every 2 seconds, until the analysis has completed or failed.
   */
  const readOutput = async (roomName: string): Promise<void> => {
    let since = 0
    let failedRequests = 0

    while (isRunning.value) {
      await wait(POLL_INTERVAL_MS)

      let response: LongPollingResponse
      try {
        response = await requestApiData<LongPollingResponse>(
          `${MANAGEMENT_LONG_POLLING_ENDPOINT}${roomName}/`,
          { method: 'GET', params: { since: String(since) } },
          MANAGEMENT_LONG_POLLING_ERROR_MESSAGE,
        )
      } catch (err: unknown) {
        failedRequests += 1
        if (failedRequests < MAX_FAILED_OUTPUT_REQUESTS) {
          continue
        }
        messages.value.push({
          level: 'error',
          text: `${MANAGEMENT_LONG_POLLING_ERROR_MESSAGE} The analysis may still be running: ${getErrorMessage(err)}`,
        })
        status.value = 'failed'
        isRunning.value = false
        return
      }

      failedRequests = 0
      messages.value.push(...response.messages)
      since = response.next
      if (response.status === 'completed' || response.status === 'failed') {
        status.value = response.status
        isRunning.value = false
      }
    }
  }

  const fetchResults = async (experimentId: number): Promise<void> => {
    const response = await requestApiData<AnalysisResultsResponse>(
      ANALYSIS_RESULTS_ENDPOINT,
      { method: 'GET', params: { experiment_id: String(experimentId) } },
      ANALYSIS_RESULTS_ERROR_MESSAGE,
    )
    results.value = response.results
  }

  /**
   * Downloads one result zip.
   *
   * Accepted data example:
   * - `experimentId = 105, name = '20260929-125825_single_Lum_CTG.zip'`
   */
  const downloadResult = async (experimentId: number, name: string): Promise<void> => {
    const zipBlob = await requestApiData<Blob>(
      ANALYSIS_DOWNLOAD_ENDPOINT,
      { method: 'GET', params: { experiment_id: String(experimentId), name }, responseType: 'blob' },
      ANALYSIS_DOWNLOAD_ERROR_MESSAGE,
    )

    const zipUrl = window.URL.createObjectURL(zipBlob)
    const link = document.createElement('a')
    link.href = zipUrl
    link.download = name
    document.body.appendChild(link)
    link.click()
    document.body.removeChild(link)
    window.URL.revokeObjectURL(zipUrl)
  }

  return {
    messages,
    status,
    isRunning,
    experimentId,
    results,
    startAnalysis,
    fetchResults,
    downloadResult,
  }
})
