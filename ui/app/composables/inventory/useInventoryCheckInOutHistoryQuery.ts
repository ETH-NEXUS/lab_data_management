import { useQuery } from '@tanstack/vue-query'
import { fetchInventoryHistory } from '~/composables/inventory/inventoryHistoryQuery.utils'
import { INVENTORY_HISTORY_QUERY_KEY } from '~/types/inventory'

export const useInventoryCheckInOutHistoryQuery = () =>
  useQuery({
    queryKey: [...INVENTORY_HISTORY_QUERY_KEY, 'check-in-out'],
    queryFn: () => fetchInventoryHistory('check_in_out'),
  })
