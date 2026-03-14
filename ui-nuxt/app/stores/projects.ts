import { defineStore } from 'pinia'
import { ref } from 'vue'
import { useAPI } from '~/composables/useAPI'
import type { PaginatedResponse } from '~/types/api'
import { PROJECTS_ENDPOINT, PROJECTS_ERROR_MESSAGE, type Project } from '~/types/projects'

const toRelativeApiEndpoint = (input: string | null | undefined): string | null => {
  if (!input) return null

  let endpoint = input
  if (endpoint.startsWith('http://') || endpoint.startsWith('https://')) {
    try {
      const parsed = new URL(endpoint)
      endpoint = `${parsed.pathname}${parsed.search}`
    } catch {
      return null
    }
  }

  return endpoint.replace(/^\/api\//, '').replace(/^\//, '')
}

export const useProjectStore = defineStore('projectStore', () => {
  const projects = ref<Project[]>([])
  const isLoading = ref(false)
  const initialized = ref(false)
  const error = ref<string | null>(null)

  const clearProjects = () => {
    projects.value = []
    error.value = null
    initialized.value = false
  }

  const loadProjects = async (force = false): Promise<Project[]> => {
    if (isLoading.value) return projects.value
    if (!force && initialized.value) return projects.value

    isLoading.value = true
    error.value = null

    const allProjects: Project[] = []
    const visitedEndpoints = new Set<string>()
    let nextEndpoint: string | null = PROJECTS_ENDPOINT

    try {
      while (nextEndpoint) {
        if (visitedEndpoints.has(nextEndpoint)) break
        visitedEndpoints.add(nextEndpoint)

        const { data, error: apiError } = await useAPI<PaginatedResponse<Project>>(nextEndpoint, {
          method: 'GET',
        })

        if (apiError.value || !data.value) {
          error.value = apiError.value?.message ?? PROJECTS_ERROR_MESSAGE
          projects.value = []
          initialized.value = true
          return projects.value
        }

        allProjects.push(...(data.value.results ?? []))
        nextEndpoint = toRelativeApiEndpoint(data.value.next)
      }

      projects.value = allProjects
      initialized.value = true
      return projects.value
    } catch (err: unknown) {
      error.value = err instanceof Error ? err.message : PROJECTS_ERROR_MESSAGE
      projects.value = []
      initialized.value = true
      return projects.value
    } finally {
      isLoading.value = false
    }
  }

  return {
    projects,
    isLoading,
    initialized,
    error,
    clearProjects,
    loadProjects,
  }
})
