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
  const { data, error } = await useAPI<ListResponse<InventoryHistoryListItem>>(
    `${INVENTORY_HISTORY_ENDPOINT}?page=1&page_size=5`,
    {
      method: 'GET',
    },
  )

  if (error.value || !data.value) {
    throw (error.value ?? new Error(INVENTORY_HISTORY_ERROR_MESSAGE)) as Error
  }

  return parseListPage(data.value).items
}

export const useInventoryHistoryQuery = () =>
  useQuery({
    queryKey: INVENTORY_HISTORY_QUERY_KEY,
    queryFn: fetchInventoryHistory,
  })
