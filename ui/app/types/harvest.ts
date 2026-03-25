export type HarvestProject = {
  id: number
  name: string
  notes?: string | null
  value?: string
  label?: string
  [key: string]: unknown
}

export type HarvestProjectsResponse = {
  projects?: HarvestProject[]
}

export type HarvestUpdateResponse = {
  success?: boolean
  [key: string]: unknown
}

export const HARVEST_PROJECTS_ENDPOINT = 'harvest/harvest_projects/'
export const getHarvestUpdateEndpoint = (projectId: number | string) => `harvest/update_harvest_info/${projectId}/`

export const HARVEST_PROJECTS_ERROR_MESSAGE = 'Failed to load Harvest projects.'
export const HARVEST_UPDATE_ERROR_MESSAGE = 'Failed to update Harvest information.'
