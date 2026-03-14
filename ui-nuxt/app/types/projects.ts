export type Project = {
  id: number
  name: string
  description?: string | null
  harvest_id?: number | null
  harvest_notes?: string | null
  created_at?: string | null
  experiments?: Array<{ id: number; name: string }>
  [key: string]: unknown
}

export const PROJECTS_QUERY_KEY = ['projects'] as const
export const PROJECTS_ENDPOINT = 'projects/'
export const PROJECTS_ERROR_MESSAGE = 'Failed to load projects.'
export const PROJECT_ERROR_MESSAGE = 'Failed to load project.'
export const PROJECT_CREATE_ERROR_MESSAGE = 'Failed to create project.'
export const PROJECT_UPDATE_ERROR_MESSAGE = 'Failed to update project.'

export const getProjectQueryKey = (projectId: number | string) => ['project', String(projectId)] as const

export type ProjectPayload = {
  name?: string
  description?: string | null
}

/**
 * Request payload for creating a project.
 *
 * Accepted payload examples:
 * - `{ name: 'Screening 2026', harvest_id: 42, harvest_notes: 'Imported from Harvest' }`
 * - `{ name: 'Manual project', harvest_id: null, harvest_notes: null }`
 */
export type CreateProjectPayload = {
  name: string
  harvest_id?: number | null
  harvest_notes?: string | null
}
