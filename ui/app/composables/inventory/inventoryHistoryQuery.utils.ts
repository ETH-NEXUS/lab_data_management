import { useAPI } from '~/composables/useAPI'
import { parseListPage, type ListResponse } from '~/composables/inventory/listResponse'
import {
  INVENTORY_HISTORY_ENDPOINT,
  INVENTORY_HISTORY_ERROR_MESSAGE,
  type InventoryHistoryListItem,
} from '~/types/inventory'

export type InventoryHistoryActivityGroup = 'check_in_out'

export const fetchInventoryHistory = async (
  activityGroup?: InventoryHistoryActivityGroup,
): Promise<InventoryHistoryListItem[]> => {
  const activityGroupQuery = activityGroup ? `&activity_group=${activityGroup}` : ''
  const { data, error } = await useAPI<ListResponse<InventoryHistoryListItem>>(
    `${INVENTORY_HISTORY_ENDPOINT}?page=1&page_size=5${activityGroupQuery}`,
    {
      method: 'GET',
    },
  )

  if (error.value || !data.value) {
    throw (error.value ?? new Error(INVENTORY_HISTORY_ERROR_MESSAGE)) as Error
  }

  return parseListPage(data.value).items
}
