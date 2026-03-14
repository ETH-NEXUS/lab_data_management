import { defineStore } from 'pinia'
import { ref } from 'vue'
import { useAPI } from '~/composables/useAPI'
import {
  PROJECT_CREATE_ERROR_MESSAGE,
  PROJECT_UPDATE_ERROR_MESSAGE,
  type CreateProjectPayload,
  type Project,
  type ProjectPayload,
} from '~/types/projects'

export const useProjectStore = defineStore('projectStore', () => {
  const isCreatingProject = ref(false)
  const isUpdatingProject = ref(false)

  /**
   * Creates a new project in the backend.
   *
   * Accepted payload examples:
   * - `{ name: 'Screening 2026', harvest_id: 10, harvest_notes: 'From Harvest' }`
   * - `{ name: 'Manual project', harvest_id: null, harvest_notes: null }`
   *
   * Returned project example:
   * - `{ id: 15, name: 'Screening 2026', experiments: [] }`
   */
  const createProject = async (payload: CreateProjectPayload): Promise<Project> => {
    isCreatingProject.value = true
    try {
      const { data, error } = await useAPI<Project>('projects/', {
        method: 'POST',
        body: payload,
      })

      if (error.value || !data.value) {
        throw (error.value ?? new Error(PROJECT_CREATE_ERROR_MESSAGE)) as Error
      }

      return data.value
    } finally {
      isCreatingProject.value = false
    }
  }

  const updateProject = async (projectId: number, payload: ProjectPayload): Promise<Project> => {
    isUpdatingProject.value = true
    try {
      const { data, error } = await useAPI<Project>(`projects/${projectId}/`, {
        method: 'PATCH',
        body: payload,
      })

      if (error.value || !data.value) {
        throw (error.value ?? new Error(PROJECT_UPDATE_ERROR_MESSAGE)) as Error
      }

      return data.value
    } finally {
      isUpdatingProject.value = false
    }
  }

  return {
    isCreatingProject,
    isUpdatingProject,
    createProject,
    updateProject,
  }
})
