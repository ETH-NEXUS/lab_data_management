import { computed, type ComputedRef } from 'vue'
import { useQuery } from '@tanstack/vue-query'
import { useAPI } from '~/composables/useAPI'
import type { PaginatedResponse } from '~/types/api'
import { toRelativeApiEndpoint } from '~/utils/apiEndpoint'
import {
  INVENTORY_STOCKS_ENDPOINT,
  INVENTORY_STOCKS_ERROR_MESSAGE,
  INVENTORY_STOCK_ERROR_MESSAGE,
  INVENTORY_STOCKS_QUERY_KEY,
  getInventoryStockQueryKey,
  type InventoryStockDetail,
  type InventoryStockListItem,
} from '~/types/inventory'

const fetchInventoryStocks = async (): Promise<InventoryStockListItem[]> => {
  const allItems: InventoryStockListItem[] = []
  const visitedEndpoints = new Set<string>()
  let nextEndpoint: string | null = INVENTORY_STOCKS_ENDPOINT

  while (nextEndpoint) {
    if (visitedEndpoints.has(nextEndpoint)) break
    visitedEndpoints.add(nextEndpoint)

    const { data, error } = await useAPI<PaginatedResponse<InventoryStockListItem>>(nextEndpoint, {
      method: 'GET',
    })

    if (error.value || !data.value) {
      throw (error.value ?? new Error(INVENTORY_STOCKS_ERROR_MESSAGE)) as Error
    }

    allItems.push(...(data.value.results ?? []))
    nextEndpoint = toRelativeApiEndpoint(data.value.next)
  }

  return allItems
}

const fetchInventoryStock = async (stockId: number): Promise<InventoryStockDetail> => {
  const { data, error } = await useAPI<InventoryStockDetail>(`inventory/stocks/${stockId}/`, {
    method: 'GET',
  })

  if (error.value || !data.value) {
    throw (error.value ?? new Error(INVENTORY_STOCK_ERROR_MESSAGE)) as Error
  }

  return data.value
}

export const useInventoryStocksQuery = () =>
  useQuery({
    queryKey: INVENTORY_STOCKS_QUERY_KEY,
    queryFn: fetchInventoryStocks,
  })

export const useInventoryStockQuery = (stockIdRef: ComputedRef<number>) =>
  useQuery({
    queryKey: computed(() => getInventoryStockQueryKey(stockIdRef.value)),
    enabled: computed(() => Number.isInteger(stockIdRef.value) && stockIdRef.value > 0),
    queryFn: () => fetchInventoryStock(stockIdRef.value),
  })
