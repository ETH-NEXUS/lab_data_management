import type { FileSystemItem } from '~/types/lab'

/**
 * API constants and payloads for management filesystem features.
 *
 * Data examples:
 * - directory response: `{ directory_content: { type: 'directory', name: 'data', path: '/data', children: [] } }`
 * - delete payload: `{ path: '/data/file.txt' }`
 * - download payload: `{ file_path: '/data/file.txt' }`
 * - file content response: `{ content: 'line 1\\nline 2' }`
 */

export const MANAGEMENT_DIRECTORY_CONTENT_ENDPOINT = 'directory_content/'
export const MANAGEMENT_RUN_COMMAND_ENDPOINT = 'run_command/'
export const MANAGEMENT_LONG_POLLING_ENDPOINT = 'long_polling/'
export const MANAGEMENT_DELETE_FILE_ENDPOINT = 'delete_file/'
export const MANAGEMENT_DOWNLOAD_FILE_ENDPOINT = 'download_file/'
export const MANAGEMENT_UPLOAD_FILE_ENDPOINT = 'upload_file/'
export const MANAGEMENT_GET_FILE_CONTENT_ENDPOINT = 'get_file_content/'

export const MANAGEMENT_DIRECTORY_CONTENT_ERROR_MESSAGE = 'Failed to load directory content.'
export const MANAGEMENT_RUN_COMMAND_ERROR_MESSAGE = 'Failed to run management command.'
export const MANAGEMENT_DELETE_FILE_ERROR_MESSAGE = 'Failed to delete file or directory.'
export const MANAGEMENT_DOWNLOAD_FILE_ERROR_MESSAGE = 'Failed to download file.'
export const MANAGEMENT_UPLOAD_FILE_ERROR_MESSAGE = 'Failed to upload file.'
export const MANAGEMENT_GET_FILE_CONTENT_ERROR_MESSAGE = 'Failed to load file content.'
export const MANAGEMENT_LONG_POLLING_ERROR_MESSAGE = 'Failed to read command output.'

export type DirectoryContentResponse = {
  directory_content: FileSystemItem
}

export type LongPollingResponse = {
  message?: string
  status?: string
}

export type FileContentResponse = {
  content: string
}
