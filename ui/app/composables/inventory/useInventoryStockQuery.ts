import { computed, type ComputedRef } from 'vue'
import { useQuery } from '@tanstack/vue-query'
import { useAPI } from '~/composables/useAPI'
import { parseListPage, type ListResponse } from '~/composables/inventory/listResponse'
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

    const { data, error } = await useAPI<ListResponse<InventoryStockListItem>>(nextEndpoint, {
      method: 'GET',
    })

    if (error.value || !data.value) {
      throw (error.value ?? new Error(INVENTORY_STOCKS_ERROR_MESSAGE)) as Error
    }

    const parsedPage = parseListPage(data.value)
    allItems.push(...parsedPage.items)
    nextEndpoint = parsedPage.nextEndpoint
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

export const useInventoryStocksQuery = (enabledRef?: ComputedRef<boolean>) =>
  useQuery({
    queryKey: INVENTORY_STOCKS_QUERY_KEY,
    enabled: computed(() => enabledRef?.value ?? true),
    queryFn: fetchInventoryStocks,
  })

export const useInventoryStockQuery = (stockIdRef: ComputedRef<number>) =>
  useQuery({
    queryKey: computed(() => getInventoryStockQueryKey(stockIdRef.value)),
    enabled: computed(() => Number.isInteger(stockIdRef.value) && stockIdRef.value > 0),
    queryFn: () => fetchInventoryStock(stockIdRef.value),
  })
