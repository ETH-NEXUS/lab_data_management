import { useQuery } from '@tanstack/vue-query'
import { useAPI } from '~/composables/useAPI'
import type { PaginatedResponse } from '~/types/api'
import { PROJECTS_ENDPOINT, PROJECTS_ERROR_MESSAGE, PROJECTS_QUERY_KEY, type Project } from '~/types/projects'

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

const fetchProjects = async (): Promise<Project[]> => {
  const allProjects: Project[] = []
  const visitedEndpoints = new Set<string>()
  let nextEndpoint: string | null = PROJECTS_ENDPOINT

  while (nextEndpoint) {
    if (visitedEndpoints.has(nextEndpoint)) break
    visitedEndpoints.add(nextEndpoint)

    const { data, error } = await useAPI<PaginatedResponse<Project>>(nextEndpoint, {
      method: 'GET',
    })

    if (error.value || !data.value) {
      throw (error.value ?? new Error(PROJECTS_ERROR_MESSAGE)) as Error
    }

    allProjects.push(...(data.value.results ?? []))
    nextEndpoint = toRelativeApiEndpoint(data.value.next)
  }

  return allProjects
}

export const useProjectsQuery = () =>
  useQuery({
    queryKey: PROJECTS_QUERY_KEY,
    queryFn: fetchProjects,
  })
