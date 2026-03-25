import { defineStore } from 'pinia'
import { ref } from 'vue'
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

export const useManagementStore = defineStore('managementStore', () => {
  const dataDirectory = ref<FileSystemItem>(createEmptyDirectoryItem())
  const selectedPath = ref('')
  const selectedPaths = ref<string[]>([])
  const commandOutput = ref('')

  const isLoadingDirectoryContent = ref(false)
  const isRunningCommand = ref(false)
  const isDeletingFile = ref(false)
  const isDownloadingFile = ref(false)
  const isUploadingFile = ref(false)
  const isLoadingFileContent = ref(false)
  const error = ref<string | null>(null)

  /**
   * Loads directory tree from backend and stores it.
   *
   * Returned data example:
   * - `{ type: 'directory', name: 'data', path: '/data', children: [{ type: 'file', name: 'x.txt', path: '/data/x.txt' }] }`
   */
  const fetchDataDirectory = async (): Promise<FileSystemItem> => {
    isLoadingDirectoryContent.value = true
    error.value = null

    try {
      const response = await requestApiData<DirectoryContentResponse>(
        MANAGEMENT_DIRECTORY_CONTENT_ENDPOINT,
        { method: 'GET' },
        MANAGEMENT_DIRECTORY_CONTENT_ERROR_MESSAGE,
      )

      dataDirectory.value = response.directory_content ?? createEmptyDirectoryItem()
      return dataDirectory.value
    } catch (err: unknown) {
      error.value = getErrorMessage(err)
      throw err
    } finally {
      isLoadingDirectoryContent.value = false
    }
  }

  /**
   * Legacy-compatible initialize action.
   */
  const initialize = async (): Promise<void> => {
    await fetchDataDirectory()
  }

  /**
   * Legacy-compatible alias used by old management page/components.
   *
   * Returned data example:
   * - `{ type: 'directory', name: 'data', path: '/data', children: [] }`
   */
  const getDataDirectory = async (): Promise<FileSystemItem> => {
    return fetchDataDirectory()
  }

  /**
   * Clears terminal-like output shown in management command cards.
   */
  const clearCommandOutput = (): void => {
    commandOutput.value = ''
  }

  /**
   * Adds one selected path to the legacy selected list.
   *
   * Accepted data example:
   * - `'/data/imports/file.csv'`
   */
  const addSelectedPath = (path: string): void => {
    selectedPath.value = path
    selectedPaths.value.push(path)
  }

  /**
   * Removes first matching selected path from the legacy selected list.
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
   * Clears all selected management paths.
   */
  const clearSelectedPaths = (): void => {
    selectedPath.value = ''
    selectedPaths.value = []
  }

  /**
   * Starts one management command and begins long-polling output.
   *
   * Accepted data example:
   * - `{ room_name: 'import_2026_03_15', command: 'echo_import' }`
   */
  const runCommand = async (formData: GeneralFormData): Promise<void> => {
    isRunningCommand.value = true
    error.value = null

    try {
      const roomNameValue = formData.room_name
      const roomName = typeof roomNameValue === 'string' ? roomNameValue : ''

      // Keep old behavior: fire polling without waiting for full completion.
      if (roomName !== '') {
        void startLongPolling(roomName)
      }

      await requestApiVoid(
        MANAGEMENT_RUN_COMMAND_ENDPOINT,
        {
          method: 'POST',
          body: { form_data: formData },
        },
        MANAGEMENT_RUN_COMMAND_ERROR_MESSAGE,
      )
    } catch (err: unknown) {
      error.value = getErrorMessage(err)
      throw err
    } finally {
      isRunningCommand.value = false
    }
  }

  /**
   * Polls command output endpoint until backend status is `completed`.
   *
   * Accepted data example:
   * - `roomName = 'import_2026_03_15'`
   */
  const startLongPolling = async (roomName: string): Promise<void> => {
    try {
      const response = await requestApiData<LongPollingResponse>(
        `${MANAGEMENT_LONG_POLLING_ENDPOINT}${roomName}/`,
        { method: 'GET' },
        MANAGEMENT_LONG_POLLING_ERROR_MESSAGE,
      )

      const message = response.message
      if (typeof message === 'string' && message !== '') {
        commandOutput.value += `\n${message}`
      }

      const status = response.status
      if (status !== 'completed') {
        setTimeout(() => {
          void startLongPolling(roomName)
        }, 300)
      } else {
        // Keep directory tree fresh after command completion.
        await fetchDataDirectory()
      }
    } catch (err: unknown) {
      // Keep old behavior: only log polling errors, do not throw to UI flow.
      console.error(MANAGEMENT_LONG_POLLING_ERROR_MESSAGE, err)
    }
  }

  /**
   * Deletes one file/directory path and refreshes directory content.
   *
   * Accepted path example:
   * - `'/data/imports/file.txt'`
   */
  const deleteFile = async (path: string): Promise<void> => {
    isDeletingFile.value = true
    error.value = null

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
    } catch (err: unknown) {
      error.value = getErrorMessage(err)
      throw err
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
    isDownloadingFile.value = true
    error.value = null

    try {
      const blob = await requestApiData<Blob>(
        MANAGEMENT_DOWNLOAD_FILE_ENDPOINT,
        {
          method: 'POST',
          body: { file_path: path },
          responseType: 'blob',
        },
        MANAGEMENT_DOWNLOAD_FILE_ERROR_MESSAGE,
      )

      return blob
    } catch (err: unknown) {
      error.value = getErrorMessage(err)
      throw err
    } finally {
      isDownloadingFile.value = false
    }
  }

  /**
   * Uploads one file into one target directory and refreshes tree.
   *
   * Accepted input examples:
   * - `directoryPath = '/data/imports', file.name = 'new.csv'`
   */
  const uploadFile = async (directoryPath: string, file: File): Promise<void> => {
    isUploadingFile.value = true
    error.value = null

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
    } catch (err: unknown) {
      error.value = getErrorMessage(err)
      throw err
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
    isLoadingFileContent.value = true
    error.value = null

    try {
      const response = await requestApiData<FileContentResponse>(
        MANAGEMENT_GET_FILE_CONTENT_ENDPOINT,
        {
          method: 'POST',
          body: { file_path: path },
        },
        MANAGEMENT_GET_FILE_CONTENT_ERROR_MESSAGE,
      )

      return response.content ?? ''
    } catch (err: unknown) {
      error.value = getErrorMessage(err)
      throw err
    } finally {
      isLoadingFileContent.value = false
    }
  }

  return {
    dataDirectory,
    selectedPath,
    selectedPaths,
    commandOutput,
    isLoadingDirectoryContent,
    isRunningCommand,
    isDeletingFile,
    isDownloadingFile,
    isUploadingFile,
    isLoadingFileContent,
    error,
    fetchDataDirectory,
    getDataDirectory,
    initialize,
    clearCommandOutput,
    addSelectedPath,
    removeSelectedPath,
    clearSelectedPaths,
    runCommand,
    startLongPolling,
    deleteFile,
    downloadFile,
    uploadFile,
    getFileContent,
  }
})
