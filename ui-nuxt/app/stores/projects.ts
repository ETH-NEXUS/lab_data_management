import { defineStore } from 'pinia'
import { ref } from 'vue'
import {
  PROJECT_CREATE_ERROR_MESSAGE,
  PROJECT_UPDATE_ERROR_MESSAGE,
  type CreateProjectPayload,
  type Project,
  type ProjectPayload,
} from '~/types/projects'
import { requestApiData } from '~/utils/apiRequests'

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
      const project = await requestApiData<Project>(
        'projects/',
        {
          method: 'POST',
          body: payload,
        },
        PROJECT_CREATE_ERROR_MESSAGE,
      )

      return project
    } finally {
      isCreatingProject.value = false
    }
  }

  const updateProject = async (projectId: number, payload: ProjectPayload): Promise<Project> => {
    isUpdatingProject.value = true
    try {
      const project = await requestApiData<Project>(
        `projects/${projectId}/`,
        {
          method: 'PATCH',
          body: payload,
        },
        PROJECT_UPDATE_ERROR_MESSAGE,
      )

      return project
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
