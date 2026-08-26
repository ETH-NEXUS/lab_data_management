import { useQuery } from '@tanstack/vue-query'
import { useAPI } from '~/composables/useAPI'
import { parseListPage, type ListResponse } from '~/composables/inventory/listResponse'
import {
  INVENTORY_HISTORY_ENDPOINT,
  INVENTORY_HISTORY_ERROR_MESSAGE,
  INVENTORY_HISTORY_QUERY_KEY,
  type InventoryHistoryListItem,
} from '~/types/inventory'

const fetchInventoryHistory = async (): Promise<InventoryHistoryListItem[]> => {
  const allItems: InventoryHistoryListItem[] = []
  const visitedEndpoints = new Set<string>()
  let nextEndpoint: string | null = INVENTORY_HISTORY_ENDPOINT

  while (nextEndpoint) {
    if (visitedEndpoints.has(nextEndpoint)) break
    visitedEndpoints.add(nextEndpoint)

    const { data, error } = await useAPI<ListResponse<InventoryHistoryListItem>>(nextEndpoint, {
      method: 'GET',
    })

    if (error.value || !data.value) {
      throw (error.value ?? new Error(INVENTORY_HISTORY_ERROR_MESSAGE)) as Error
    }

    const parsedPage = parseListPage(data.value)
    allItems.push(...parsedPage.items)
    nextEndpoint = parsedPage.nextEndpoint
  }

  return allItems
}

export const useInventoryHistoryQuery = () =>
  useQuery({
    queryKey: INVENTORY_HISTORY_QUERY_KEY,
    queryFn: fetchInventoryHistory,
  })
