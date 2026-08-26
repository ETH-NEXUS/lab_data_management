import { useAPI } from '~/composables/useAPI'
import {
  INVENTORY_HISTORY_ENDPOINT,
  INVENTORY_HISTORY_ERROR_MESSAGE,
  type InventoryHistoryListItem,
} from '~/types/inventory'
import type { PaginatedResponse } from '~/types/api'

export type InventoryHistoryActivityGroup = 'check_in_out'

export type InventoryHistoryPageQueryParams = {
  page: number
  pageSize: number
  activityGroup?: InventoryHistoryActivityGroup
}

/**
 * Loads one server-side page of inventory history.
 *
 * Request example:
 * - `{ page: 2, pageSize: 20 }`
 * - `{ page: 1, pageSize: 5, activityGroup: 'check_in_out' }`
 */
export const fetchInventoryHistoryPage = async (
  params: InventoryHistoryPageQueryParams,
): Promise<PaginatedResponse<InventoryHistoryListItem>> => {
  const requestParams: Record<string, string> = {
    page: String(params.page),
    page_size: String(params.pageSize),
  }

  if (params.activityGroup) {
    requestParams.activity_group = params.activityGroup
  }

  const { data, error } = await useAPI<PaginatedResponse<InventoryHistoryListItem>>(INVENTORY_HISTORY_ENDPOINT, {
    method: 'GET',
    params: requestParams,
  })

  if (error.value || !data.value) {
    throw (error.value ?? new Error(INVENTORY_HISTORY_ERROR_MESSAGE)) as Error
  }

  return data.value
}

export const fetchInventoryHistory = async (
  activityGroup?: InventoryHistoryActivityGroup,
): Promise<InventoryHistoryListItem[]> => {
  const historyPage = await fetchInventoryHistoryPage({
    page: 1,
    pageSize: 5,
    activityGroup,
  })

  return historyPage.results ?? []
}
