import { defineStore } from 'pinia'
import { computed, ref } from 'vue'
import type { FileSystemItem, GeneralFormData } from '~/types/lab'
import {
  MANAGEMENT_DELETE_FILE_ENDPOINT,
  MANAGEMENT_DELETE_FILE_ERROR_MESSAGE,
  MANAGEMENT_DIRECTORY_CONTENT_ENDPOINT,
  MANAGEMENT_DIRECTORY_CONTENT_ERROR_MESSAGE,
  MANAGEMENT_DOWNLOAD_FILE_ENDPOINT,
  MANAGEMENT_DOWNLOAD_FILE_ERROR_MESSAGE,
  MANAGEMENT_GET_FILE_CONTENT_ENDPOINT,
  MANAGEMENT_GET_FILE_CONTENT_ERROR_MESSAGE,
  MANAGEMENT_LONG_POLLING_ENDPOINT,
  MANAGEMENT_LONG_POLLING_ERROR_MESSAGE,
  MANAGEMENT_RUN_COMMAND_ENDPOINT,
  MANAGEMENT_RUN_COMMAND_ERROR_MESSAGE,
  MANAGEMENT_UPLOAD_FILE_ENDPOINT,
  MANAGEMENT_UPLOAD_FILE_ERROR_MESSAGE,
  type CommandMessage,
  type CommandStatus,
  type DirectoryContentResponse,
  type FileContentResponse,
  type LongPollingResponse,
} from '~/types/management'
import { requestApiData, requestApiVoid } from '~/utils/apiRequests'
import { getErrorMessage } from '~/utils/errors'

const createEmptyDirectoryItem = (): FileSystemItem => ({
  type: 'directory',
  name: '',
  path: '',
  children: [],
})

// The command output is asked for this many times in a row before giving up
const MAX_FAILED_OUTPUT_REQUESTS = 5
const FAILED_OUTPUT_REQUEST_DELAY_MS = 1000
// A command that has not written a single line in this time has not started yet:
// it waits in the queue behind other commands, or the worker is down
const COMMAND_START_TIMEOUT_MS = 60000
// Shown once, as info: a command that starts later still ends as usual, and a
// warning would mark it "completed with warnings" although nothing was wrong
const COMMAND_WAITING_MESSAGE =
  'The command has not started yet: it is waiting for the worker (other commands run first), or the worker is down.'
// How often the page asks for new output while a command runs
const OUTPUT_REQUEST_DELAY_MS = 1000
// The command of this browser, so its output comes back after a reload
const ROOM_NAME_STORAGE_KEY = 'managementRoomName'

export const useManagementStore = defineStore('managementStore', () => {
  const dataDirectory = ref<FileSystemItem>(createEmptyDirectoryItem())
  const selectedPath = ref('')
  const selectedPaths = ref<string[]>([])
  const commandMessages = ref<CommandMessage[]>([])
  const commandStatus = ref<CommandStatus | null>(null)
  // Polling reads only the output of this room; a new command replaces it
  const activeRoomName = ref('')
  // Every start of the output raises this number, so an older loop stops
  const pollRound = ref(0)
  const commandRequestError = ref<string | null>(null)

  const isLoadingDirectoryContent = ref(false)
  const isRunningCommand = ref(false)
  const isDeletingFile = ref(false)
  const isUploadingFile = ref(false)

  /**
   * Loads directory tree from backend and stores it.
   *
   * Returned data example:
   * - `{ type: 'directory', name: 'data', path: '/data', children: [{ type: 'file', name: 'x.txt', path: '/data/x.txt' }] }`
   */
  const fetchDataDirectory = async (): Promise<FileSystemItem> => {
    isLoadingDirectoryContent.value = true

    try {
      const response = await requestApiData<DirectoryContentResponse>(
        MANAGEMENT_DIRECTORY_CONTENT_ENDPOINT,
        { method: 'GET' },
        MANAGEMENT_DIRECTORY_CONTENT_ERROR_MESSAGE,
      )

      dataDirectory.value = response.directory_content ?? createEmptyDirectoryItem()
      return dataDirectory.value
    } finally {
      isLoadingDirectoryContent.value = false
    }
  }

  /**
   * Loads the directory tree when the page opens. The page and the navigation
   * tree both call this, but the folders are only read once.
   */
  const initialize = async (): Promise<void> => {
    if (dataDirectory.value.children.length > 0) {
      return
    }
    await fetchDataDirectory()
  }

  /**
   * The command this browser started last, so its output comes back after a
   * page reload. Browser storage can be switched off, then nothing is remembered.
   */
  const rememberedRoomName = computed({
    get: (): string => {
      try {
        return localStorage.getItem(ROOM_NAME_STORAGE_KEY) ?? ''
      } catch {
        return ''
      }
    },
    set: (roomName: string): void => {
      try {
        if (roomName === '') {
          localStorage.removeItem(ROOM_NAME_STORAGE_KEY)
        } else {
          localStorage.setItem(ROOM_NAME_STORAGE_KEY, roomName)
        }
      } catch {
        // Nothing is remembered, the output is only lost after a reload
      }
    },
  })

  /**
   * Shows the output of the command of this browser again, when the page is
   * opened while it runs or after it has ended. Called once per page.
   */
  const resumeCommandOutput = async (): Promise<void> => {
    const roomName = rememberedRoomName.value
    // Nothing was started from this browser, or this page runs a command itself
    if (roomName === '' || activeRoomName.value !== '') {
      return
    }
    const round = pollRound.value

    let response: LongPollingResponse
    try {
      response = await requestApiData<LongPollingResponse>(
        `${MANAGEMENT_LONG_POLLING_ENDPOINT}${roomName}/`,
        { method: 'GET', params: { since: '0' } },
        MANAGEMENT_LONG_POLLING_ERROR_MESSAGE,
      )
    } catch (err: unknown) {
      // The command stays remembered, so opening the page again tries once more
      console.error(MANAGEMENT_LONG_POLLING_ERROR_MESSAGE, err)
      return
    }

    // A command that was started while the answer was on its way keeps the page
    if (round !== pollRound.value || activeRoomName.value !== '') {
      return
    }

    commandMessages.value = response.messages
    commandStatus.value = response.status
    // The output of a command that has ended is shown once
    rememberedRoomName.value = ''

    if (response.status !== 'running') {
      return
    }

    activeRoomName.value = roomName
    rememberedRoomName.value = roomName
    isRunningCommand.value = true
    void pollCommandOutput(roomName, response.next, 0, Date.now(), round)
  }

  /**
   * Removes one path from the list of selected paths.
   *
   * Accepted data example:
   * - `'/data/imports/file.csv'`
   */
  const removeSelectedPath = (path: string): void => {
    const index = selectedPaths.value.indexOf(path)
    if (index >= 0) {
      selectedPaths.value.splice(index, 1)
    }
  }

  /**
   * Starts one management command and reads its output while it runs.
   * The request only starts the command (it runs in the celery container), so
   * `isRunningCommand` stays true until the output says the command has ended.
   *
   * Accepted data example:
   * - `{ room_name: '12_1726563600000', command: 'map', machine: 'echo', path: '/data/run_1' }`
   */
  const runCommand = async (formData: GeneralFormData): Promise<void> => {
    isRunningCommand.value = true
    // A loop that still reads the output of an earlier command stops now
    pollRound.value += 1

    const roomName = typeof formData.room_name === 'string' ? formData.room_name : ''
    commandMessages.value = [{ level: 'info', text: `Executing command: ${String(formData.command)}` }]
    commandStatus.value = 'running'
    commandRequestError.value = null
    activeRoomName.value = roomName
    rememberedRoomName.value = roomName

    if (roomName !== '') {
      void pollCommandOutput(roomName, 0, 0, Date.now(), pollRound.value)
    }

    try {
      await requestApiVoid(
        MANAGEMENT_RUN_COMMAND_ENDPOINT,
        {
          method: 'POST',
          body: { form_data: formData },
        },
        MANAGEMENT_RUN_COMMAND_ERROR_MESSAGE,
      )
    } catch (err: unknown) {
      // The command did not start, polling stops after its next read
      commandRequestError.value = getErrorMessage(err)
      isRunningCommand.value = false
    }

    // Without a room name there is no output to wait for
    if (roomName === '') {
      isRunningCommand.value = false
    }
  }

  /**
   * Adds the new output of a command every second, until the command has
   * completed or failed. A failed request is repeated, up to 5 requests in a row.
   * A loop of an earlier command (`round`) stops as soon as a new one starts.
   *
   * Accepted data example:
   * - `roomName = '12_1726563600000', since = 4` (the first 4 messages are already shown)
   */
  const pollCommandOutput = async (
    roomName: string,
    since: number,
    failedRequests = 0,
    startedAt = Date.now(),
    round = pollRound.value,
  ): Promise<void> => {
    if (round !== pollRound.value || roomName !== activeRoomName.value) {
      return
    }

    let response: LongPollingResponse
    try {
      response = await requestApiData<LongPollingResponse>(
        `${MANAGEMENT_LONG_POLLING_ENDPOINT}${roomName}/`,
        { method: 'GET', params: { since: String(since) } },
        MANAGEMENT_LONG_POLLING_ERROR_MESSAGE,
      )
    } catch (err: unknown) {
      console.error(MANAGEMENT_LONG_POLLING_ERROR_MESSAGE, err)
      // A short network problem must not stop showing the output, so ask again
      if (failedRequests + 1 < MAX_FAILED_OUTPUT_REQUESTS) {
        setTimeout(() => {
          void pollCommandOutput(roomName, since, failedRequests + 1, startedAt, round)
        }, FAILED_OUTPUT_REQUEST_DELAY_MS)
        return
      }
      commandMessages.value.push({
        level: 'error',
        text: `${MANAGEMENT_LONG_POLLING_ERROR_MESSAGE} The command may still be running.`,
      })
      commandStatus.value = 'failed'
      isRunningCommand.value = false
      return
    }

    if (round !== pollRound.value || roomName !== activeRoomName.value) {
      return
    }
    commandMessages.value.push(...response.messages)

    if (commandRequestError.value !== null) {
      commandMessages.value.push({ level: 'error', text: commandRequestError.value })
      commandStatus.value = 'failed'
      return
    }

    if (response.status === 'completed' || response.status === 'failed') {
      commandStatus.value = response.status
      isRunningCommand.value = false
      rememberedRoomName.value = ''
      // The command may have created or changed files
      try {
        await fetchDataDirectory()
      } catch (err: unknown) {
        console.error(MANAGEMENT_DIRECTORY_CONTENT_ERROR_MESSAGE, err)
      }
      return
    }

    // Without a single line after a minute the command waits for the worker. The
    // page keeps asking, so the output shows up as soon as the worker takes it.
    const hasNotStarted = response.next === 0 && Date.now() - startedAt > COMMAND_START_TIMEOUT_MS
    const waitingIsShown = commandMessages.value.some((message) => message.text === COMMAND_WAITING_MESSAGE)
    if (hasNotStarted && !waitingIsShown) {
      commandMessages.value.push({ level: 'info', text: COMMAND_WAITING_MESSAGE })
    }

    setTimeout(() => {
      void pollCommandOutput(roomName, response.next, 0, startedAt, round)
    }, OUTPUT_REQUEST_DELAY_MS)
  }

  /**
   * Deletes one file/directory path and refreshes directory content.
   *
   * Accepted path example:
   * - `'/data/imports/file.txt'`
   */
  const deleteFile = async (path: string): Promise<void> => {
    isDeletingFile.value = true

    try {
      await requestApiVoid(
        MANAGEMENT_DELETE_FILE_ENDPOINT,
        {
          method: 'POST',
          body: { path },
        },
        MANAGEMENT_DELETE_FILE_ERROR_MESSAGE,
      )

      await fetchDataDirectory()
    } finally {
      isDeletingFile.value = false
    }
  }

  /**
   * Downloads one file as Blob.
   *
   * Accepted path example:
   * - `'/data/imports/file.txt'`
   */
  const downloadFile = async (path: string): Promise<Blob> => {
    return await requestApiData<Blob>(
      MANAGEMENT_DOWNLOAD_FILE_ENDPOINT,
      {
        method: 'POST',
        body: { file_path: path },
        responseType: 'blob',
      },
      MANAGEMENT_DOWNLOAD_FILE_ERROR_MESSAGE,
    )
  }

  /**
   * Uploads one file into one target directory and refreshes tree.
   *
   * Accepted input examples:
   * - `directoryPath = '/data/imports', file.name = 'new.csv'`
   */
  const uploadFile = async (directoryPath: string, file: File): Promise<void> => {
    isUploadingFile.value = true

    try {
      const formData = new FormData()
      formData.append('file', file)
      formData.append('directory_path', directoryPath)

      await requestApiVoid(
        MANAGEMENT_UPLOAD_FILE_ENDPOINT,
        {
          method: 'POST',
          body: formData,
        },
        MANAGEMENT_UPLOAD_FILE_ERROR_MESSAGE,
      )

      await fetchDataDirectory()
    } finally {
      isUploadingFile.value = false
    }
  }

  /**
   * Reads text content of one file path.
   *
   * Accepted path example:
   * - `'/data/imports/readme.txt'`
   *
   * Returned data example:
   * - `'first line\\nsecond line'`
   */
  const getFileContent = async (path: string): Promise<string> => {
    const response = await requestApiData<FileContentResponse>(
      MANAGEMENT_GET_FILE_CONTENT_ENDPOINT,
      {
        method: 'POST',
        body: { file_path: path },
      },
      MANAGEMENT_GET_FILE_CONTENT_ERROR_MESSAGE,
    )

    return response.content ?? ''
  }

  return {
    dataDirectory,
    selectedPath,
    selectedPaths,
    commandMessages,
    commandStatus,
    isLoadingDirectoryContent,
    isRunningCommand,
    isDeletingFile,
    isUploadingFile,
    fetchDataDirectory,
    initialize,
    resumeCommandOutput,
    removeSelectedPath,
    runCommand,
    deleteFile,
    downloadFile,
    uploadFile,
    getFileContent,
  }
})
