import { computed, type ComputedRef } from 'vue'
import { useQuery } from '@tanstack/vue-query'
import { useAPI } from '~/composables/useAPI'
import { parseListPage, type ListResponse } from '~/composables/inventory/listResponse'
import {
  INVENTORY_USAGES_ENDPOINT,
  INVENTORY_USAGES_ERROR_MESSAGE,
  INVENTORY_USAGE_ERROR_MESSAGE,
  INVENTORY_USAGES_QUERY_KEY,
  getInventoryUsageQueryKey,
  type InventoryUsageDetail,
  type InventoryUsageListItem,
} from '~/types/inventory'

const fetchInventoryUsages = async (): Promise<InventoryUsageListItem[]> => {
  const allItems: InventoryUsageListItem[] = []
  const visitedEndpoints = new Set<string>()
  let nextEndpoint: string | null = INVENTORY_USAGES_ENDPOINT

  while (nextEndpoint) {
    if (visitedEndpoints.has(nextEndpoint)) break
    visitedEndpoints.add(nextEndpoint)

    const { data, error } = await useAPI<ListResponse<InventoryUsageListItem>>(nextEndpoint, {
      method: 'GET',
    })

    if (error.value || !data.value) {
      throw (error.value ?? new Error(INVENTORY_USAGES_ERROR_MESSAGE)) as Error
    }

    const parsedPage = parseListPage(data.value)
    allItems.push(...parsedPage.items)
    nextEndpoint = parsedPage.nextEndpoint
  }

  return allItems
}

const fetchInventoryUsage = async (usageId: number): Promise<InventoryUsageDetail> => {
  const { data, error } = await useAPI<InventoryUsageDetail>(`inventory/material-usages/${usageId}/`, {
    method: 'GET',
  })

  if (error.value || !data.value) {
    throw (error.value ?? new Error(INVENTORY_USAGE_ERROR_MESSAGE)) as Error
  }

  return data.value
}

const fetchRecentProjectUsages = async (): Promise<InventoryUsageListItem[]> => {
  const { data, error } = await useAPI<InventoryUsageListItem[]>(`${INVENTORY_USAGES_ENDPOINT}recent_project/`, {
    method: 'GET',
  })

  if (error.value || !data.value) {
    throw (error.value ?? new Error(INVENTORY_USAGES_ERROR_MESSAGE)) as Error
  }

  return data.value
}

const fetchRecentExperimentUsages = async (): Promise<InventoryUsageListItem[]> => {
  const { data, error } = await useAPI<InventoryUsageListItem[]>(`${INVENTORY_USAGES_ENDPOINT}recent_experiment/`, {
    method: 'GET',
  })

  if (error.value || !data.value) {
    throw (error.value ?? new Error(INVENTORY_USAGES_ERROR_MESSAGE)) as Error
  }

  return data.value
}

export const useInventoryUsagesQuery = () =>
  useQuery({
    queryKey: INVENTORY_USAGES_QUERY_KEY,
    queryFn: fetchInventoryUsages,
  })

export const useInventoryRecentProjectUsagesQuery = (enabledRef?: ComputedRef<boolean>) =>
  useQuery({
    queryKey: [...INVENTORY_USAGES_QUERY_KEY, 'recent-project'],
    enabled: computed(() => enabledRef?.value ?? true),
    queryFn: fetchRecentProjectUsages,
  })

export const useInventoryRecentExperimentUsagesQuery = (enabledRef?: ComputedRef<boolean>) =>
  useQuery({
    queryKey: [...INVENTORY_USAGES_QUERY_KEY, 'recent-experiment'],
    enabled: computed(() => enabledRef?.value ?? true),
    queryFn: fetchRecentExperimentUsages,
  })

export const useInventoryUsageQuery = (usageIdRef: ComputedRef<number>) =>
  useQuery({
    queryKey: computed(() => getInventoryUsageQueryKey(usageIdRef.value)),
    enabled: computed(() => Number.isInteger(usageIdRef.value) && usageIdRef.value > 0),
    queryFn: () => fetchInventoryUsage(usageIdRef.value),
  })
