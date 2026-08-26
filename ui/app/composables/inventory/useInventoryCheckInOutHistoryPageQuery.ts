import { computed, type ComputedRef } from 'vue'
import { useInventoryHistoryPageQuery } from '~/composables/inventory/useInventoryHistoryPageQuery'

export type InventoryCheckInOutHistoryPageQueryParams = {
  page: number
  pageSize: number
}

/**
 * Loads one history page containing only inventory check-ins and check-outs.
 *
 * Params example: `{ page: 1, pageSize: 20 }`.
 */
export const useInventoryCheckInOutHistoryPageQuery = (
  paramsRef: ComputedRef<InventoryCheckInOutHistoryPageQueryParams>,
) => {
  const historyPageQueryParams = computed(() => ({
    ...paramsRef.value,
    activityGroup: 'check_in_out' as const,
  }))

  return useInventoryHistoryPageQuery(historyPageQueryParams)
}
