import { computed, type ComputedRef } from 'vue'
import { useQuery } from '@tanstack/vue-query'
import {
  fetchInventoryHistoryPage,
  type InventoryHistoryPageQueryParams,
} from '~/composables/inventory/inventoryHistoryQuery.utils'
import { INVENTORY_HISTORY_QUERY_KEY } from '~/types/inventory'

/**
 * Keeps paginated history data in the query cache per requested page.
 *
 * Params example: `{ page: 1, pageSize: 20 }`.
 */
export const useInventoryHistoryPageQuery = (paramsRef: ComputedRef<InventoryHistoryPageQueryParams>) =>
  useQuery({
    queryKey: computed(() => [...INVENTORY_HISTORY_QUERY_KEY, 'page', paramsRef.value]),
    queryFn: () => fetchInventoryHistoryPage(paramsRef.value),
  })
