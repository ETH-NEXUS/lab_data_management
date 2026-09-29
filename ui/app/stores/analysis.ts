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
// An analysis without a single message after this time has not started yet.
// It still runs when the worker takes it, so the page keeps waiting (as for a
// command of the management page) and only says why nothing happens.
const START_TIMEOUT_MS = 60000
const WAITING_MESSAGE =
  'The analysis has not started yet: the analysis worker (container celery-analysis) runs one analysis at a time and another one runs first, or the worker is down.'

// The run this browser started last, so its output comes back after a page reload
const ROOM_NAME_STORAGE_KEY = 'analysis_room_name'

const wait = (milliseconds: number) => new Promise((resolve) => setTimeout(resolve, milliseconds))

const rememberRoomName = (roomName: string): void => {
  try {
    localStorage.setItem(ROOM_NAME_STORAGE_KEY, roomName)
  } catch {
    // Browser storage is switched off: the output is only lost after a reload
  }
}

const rememberedRoomName = (): string => {
  try {
    return localStorage.getItem(ROOM_NAME_STORAGE_KEY) ?? ''
  } catch {
    return ''
  }
}

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
  // The experiment whose results were asked for last
  let resultsExperimentId: number | null = null

  /**
   * Starts one analysis and shows its messages until it has ended.
   *
   * Accepted data example:
   * - `{ experiment_id: 105, label: 'Lum_CTG', analysis_type: 'single', settings: {} }`
   */
  const startAnalysis = async (payload: Omit<StartAnalysisPayload, 'room_name'>): Promise<void> => {
    // A new room for every run, so the output of an earlier run is never shown
    const roomName = `${payload.experiment_id}_${Date.now()}`
    rememberRoomName(roomName)
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
   * Shows the last analysis of this browser again after a page reload: its
   * messages so far and, while it still runs, the new ones. The room name starts
   * with the experiment id, e.g. '105_1727600000000'.
   */
  const resumeAnalysis = async (): Promise<void> => {
    const roomName = rememberedRoomName()
    // Nothing was started from this browser, or this page shows a run already
    if (roomName === '' || status.value !== null) return

    let response: LongPollingResponse
    try {
      response = await requestApiData<LongPollingResponse>(
        `${MANAGEMENT_LONG_POLLING_ENDPOINT}${roomName}/`,
        { method: 'GET', params: { since: '0' } },
        MANAGEMENT_LONG_POLLING_ERROR_MESSAGE,
      )
    } catch {
      return
    }
    // A run started while the answer was on its way keeps the page. The output
    // of a run is kept for a day; after that the status is empty.
    if (status.value !== null || response.status === null) return

    experimentId.value = Number(roomName.split('_')[0])
    messages.value = response.messages
    status.value = response.status
    if (response.status === 'running') {
      isRunning.value = true
      await readOutput(roomName, response.next)
    }
  }

  /**
   * Adds the new messages every 2 seconds, until the analysis has completed or failed.
   * `since` is the number of messages that are shown already.
   */
  const readOutput = async (roomName: string, since = 0): Promise<void> => {
    let failedRequests = 0
    const startedAt = Date.now()

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

      const hasNotStarted = since === 0 && Date.now() - startedAt > START_TIMEOUT_MS
      const waitingIsShown = messages.value.some((message) => message.text === WAITING_MESSAGE)
      if (hasNotStarted && !waitingIsShown) {
        messages.value.push({ level: 'info', text: WAITING_MESSAGE })
      }
      if (response.status === 'completed' || response.status === 'failed') {
        status.value = response.status
        isRunning.value = false
      }
    }
  }

  /**
   * Loads the result zips of one experiment. When another experiment is opened
   * before the answer arrives, the late answer is ignored, so a page never
   * lists the results of another experiment.
   */
  const fetchResults = async (experimentId: number): Promise<void> => {
    resultsExperimentId = experimentId
    const response = await requestApiData<AnalysisResultsResponse>(
      ANALYSIS_RESULTS_ENDPOINT,
      { method: 'GET', params: { experiment_id: String(experimentId) } },
      ANALYSIS_RESULTS_ERROR_MESSAGE,
    )
    if (experimentId !== resultsExperimentId) return
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
    resumeAnalysis,
    fetchResults,
    downloadResult,
  }
})
