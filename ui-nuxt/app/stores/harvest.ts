import { defineStore } from 'pinia'
import { ref } from 'vue'
import { useStorage } from '@vueuse/core'
import { useAPI } from '~/composables/useAPI'
import {
  getHarvestUpdateEndpoint,
  HARVEST_PROJECTS_ENDPOINT,
  HARVEST_PROJECTS_ERROR_MESSAGE,
  HARVEST_UPDATE_ERROR_MESSAGE,
  type HarvestProject,
  type HarvestProjectsResponse,
  type HarvestUpdateResponse,
} from '~/types/harvest'

const sortHarvestProjects = (projects: HarvestProject[]): HarvestProject[] => {
  const currentYear = new Date().getFullYear().toString()

  return [...projects].sort((a, b) => {
    const aContainsYear = a.name.includes(currentYear)
    const bContainsYear = b.name.includes(currentYear)

    if (aContainsYear && !bContainsYear) return -1
    if (!aContainsYear && bContainsYear) return 1
    return a.name.localeCompare(b.name)
  })
}

const mapHarvestProjects = (projects: HarvestProject[]): HarvestProject[] =>
  projects.map((project) => ({
    ...project,
    value: project.name,
    label: project.name,
  }))

export const useHarvestStore = defineStore('harvestStore', () => {
  // Persist list across reloads (replacement for pinia-persist in this Nuxt setup).
  const harvestProjects = useStorage<HarvestProject[]>('harvest-projects', [])

  const isLoading = ref(false)
  const isUpdatingInfo = ref(false)
  const error = ref<string | null>(null)

  const getHarvestProjects = async (): Promise<HarvestProject[]> => {
    isLoading.value = true
    error.value = null

    try {
      const { data, error: apiError } = await useAPI<HarvestProjectsResponse>(HARVEST_PROJECTS_ENDPOINT, {
        method: 'GET',
      })

      if (apiError.value || !data.value) {
        throw (apiError.value ?? new Error(HARVEST_PROJECTS_ERROR_MESSAGE)) as Error
      }

      const preparedProjects = mapHarvestProjects(data.value.projects ?? [])
      harvestProjects.value = sortHarvestProjects(preparedProjects)
      return harvestProjects.value
    } catch (err: unknown) {
      error.value = err instanceof Error ? err.message : HARVEST_PROJECTS_ERROR_MESSAGE
      throw err
    } finally {
      isLoading.value = false
    }
  }

  const updateHarvestInfo = async (projectId: number): Promise<HarvestUpdateResponse> => {
    isUpdatingInfo.value = true

    try {
      const { data, error: apiError } = await useAPI<HarvestUpdateResponse>(getHarvestUpdateEndpoint(projectId), {
        method: 'GET',
      })

      if (apiError.value || !data.value) {
        throw (apiError.value ?? new Error(HARVEST_UPDATE_ERROR_MESSAGE)) as Error
      }

      await initialize()
      return data.value
    } finally {
      isUpdatingInfo.value = false
    }
  }

  const initialize = async (): Promise<void> => {
    if (harvestProjects.value.length === 0) {
      await getHarvestProjects()
    }
  }

  return {
    harvestProjects,
    isLoading,
    isUpdatingInfo,
    error,
    getHarvestProjects,
    updateHarvestInfo,
    initialize,
  }
})
