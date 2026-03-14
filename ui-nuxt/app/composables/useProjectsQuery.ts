import { computed, type ComputedRef } from 'vue'
import { useQuery } from '@tanstack/vue-query'
import { useAPI } from '~/composables/useAPI'
import type { PaginatedResponse } from '~/types/api'
import { toRelativeApiEndpoint } from '~/utils/apiEndpoint'
import {
  getProjectQueryKey,
  PROJECT_ERROR_MESSAGE,
  PROJECTS_ENDPOINT,
  PROJECTS_ERROR_MESSAGE,
  PROJECTS_QUERY_KEY,
  type Project,
} from '~/types/projects'

/**
 * Fetches all project pages and returns a single flattened list.
 *
 * We use manual pagination here (instead of one request) because the endpoint is paginated.
 * Vue Query wraps this function and caches the final aggregated array.
 *
 * Returned list example:
 * `[{ id: 1, name: 'Nexus' }, { id: 2, name: 'Screening 2026' }]`
 */
const fetchProjects = async (): Promise<Project[]> => {
  const allProjects: Project[] = []
  // Safety guard: avoids infinite loops if backend returns a repeated `next` URL.
  const visitedEndpoints = new Set<string>()
  let nextEndpoint: string | null = PROJECTS_ENDPOINT

  while (nextEndpoint) {
    if (visitedEndpoints.has(nextEndpoint)) break
    visitedEndpoints.add(nextEndpoint)

    // Request one page at a time until backend returns `next = null`.
    const { data, error } = await useAPI<PaginatedResponse<Project>>(nextEndpoint, {
      method: 'GET',
    })

    // Vue Query expects thrown errors to set `isError` + `error`.
    if (error.value || !data.value) {
      throw (error.value ?? new Error(PROJECTS_ERROR_MESSAGE)) as Error
    }

    // Merge current page into one list consumed by navigation/UI.
    allProjects.push(...(data.value.results ?? []))
    // Normalize next-page pointer for the next loop iteration.
    nextEndpoint = toRelativeApiEndpoint(data.value.next)
  }

  return allProjects
}

/**
 * Loads one project from the API by id.
 *
 * Accepted input example:
 * - `projectId = 12`
 *
 * Returned data example:
 * - `{ id: 12, name: 'Screening 2026', description: 'Primary assay run', experiments: [{ id: 5, name: 'Exp A' }] }`
 */
const fetchProject = async (projectId: number): Promise<Project> => {
  const { data, error } = await useAPI<Project>(`projects/${projectId}/`, {
    method: 'GET',
  })

  if (error.value || !data.value) {
    throw (error.value ?? new Error(PROJECT_ERROR_MESSAGE)) as Error
  }

  return data.value
}

/**
 * Shared projects query hook.
 *
 * The stable query key allows:
 * - cache reuse across components,
 * - targeted invalidation (e.g. from refresh button),
 * - automatic status handling (loading/error/success).
 */
export const useProjectsQuery = () =>
  useQuery({
    queryKey: PROJECTS_QUERY_KEY,
    queryFn: fetchProjects,
  })

/**
 * Shared query for one project.
 *
 * Accepted input example:
 * - `projectIdRef.value = 12`
 *
 * Query key example:
 * - `['project', '12']`
 */
export const useProjectQuery = (projectIdRef: ComputedRef<number>) =>
  useQuery({
    queryKey: computed(() => getProjectQueryKey(projectIdRef.value)),
    enabled: computed(() => Number.isInteger(projectIdRef.value) && projectIdRef.value > 0),
    queryFn: () => fetchProject(projectIdRef.value),
  })
