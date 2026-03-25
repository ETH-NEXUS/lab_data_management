import { computed, type ComputedRef } from 'vue'
import { useQuery } from '@tanstack/vue-query'
import { useAPI } from '~/composables/useAPI'
import type { PaginatedResponse } from '~/types/api'
import { toRelativeApiEndpoint } from '~/utils/apiEndpoint'
import {
  EXPERIMENT_ERROR_MESSAGE,
  EXPERIMENTS_ENDPOINT,
  EXPERIMENTS_ERROR_MESSAGE,
  EXPERIMENTS_QUERY_KEY,
  getExperimentQueryKey,
  type Experiment,
} from '~/types/experiments'

/**
 * Fetches all experiment pages and returns one flattened list.
 *
 * Returned data example:
 * - `[{ id: 21, name: 'Dose response', project: 5 }, { id: 22, name: 'QC run', project: 5 }]`
 */
const fetchExperiments = async (): Promise<Experiment[]> => {
  const allExperiments: Experiment[] = []
  const visitedEndpoints = new Set<string>()
  let nextEndpoint: string | null = EXPERIMENTS_ENDPOINT

  while (nextEndpoint) {
    if (visitedEndpoints.has(nextEndpoint)) break
    visitedEndpoints.add(nextEndpoint)

    const { data, error } = await useAPI<PaginatedResponse<Experiment>>(nextEndpoint, {
      method: 'GET',
    })

    if (error.value || !data.value) {
      throw (error.value ?? new Error(EXPERIMENTS_ERROR_MESSAGE)) as Error
    }

    allExperiments.push(...(data.value.results ?? []))
    nextEndpoint = toRelativeApiEndpoint(data.value.next)
  }

  return allExperiments
}

/**
 * Fetches one experiment by id.
 *
 * Accepted input example:
 * - `experimentId = 21`
 *
 * Returned data example:
 * - `{ id: 21, name: 'Dose response', project: 5 }`
 */
const fetchExperiment = async (experimentId: number): Promise<Experiment> => {
  const { data, error } = await useAPI<Experiment>(`experiments/${experimentId}/`, {
    method: 'GET',
  })

  if (error.value || !data.value) {
    throw (error.value ?? new Error(EXPERIMENT_ERROR_MESSAGE)) as Error
  }

  return data.value
}

/**
 * Shared Vue Query hook for the full experiments list.
 */
export const useExperimentsQuery = () =>
  useQuery({
    queryKey: EXPERIMENTS_QUERY_KEY,
    queryFn: fetchExperiments,
  })

/**
 * Shared Vue Query hook for one experiment.
 *
 * Accepted input example:
 * - `experimentIdRef.value = 21`
 */
export const useExperimentQuery = (experimentIdRef: ComputedRef<number>) =>
  useQuery({
    queryKey: computed(() => getExperimentQueryKey(experimentIdRef.value)),
    enabled: computed(() => Number.isInteger(experimentIdRef.value) && experimentIdRef.value > 0),
    queryFn: () => fetchExperiment(experimentIdRef.value),
  })
