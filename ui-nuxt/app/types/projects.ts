export type Project = {
  id: number
  name: string
  harvest_id?: number | null
  [key: string]: unknown
}

export const PROJECTS_QUERY_KEY = ['projects'] as const
export const PROJECTS_ENDPOINT = 'projects/'
export const PROJECTS_ERROR_MESSAGE = 'Failed to load projects.'
