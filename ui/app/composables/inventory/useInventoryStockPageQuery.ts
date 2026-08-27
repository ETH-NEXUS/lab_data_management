import { computed, type ComputedRef } from 'vue'
import { useQuery } from '@tanstack/vue-query'
import { useAPI } from '~/composables/useAPI'
import {
  INVENTORY_STOCKS_ENDPOINT,
  INVENTORY_STOCKS_ERROR_MESSAGE,
  type InventoryStockListItem,
  type InventoryStockPreset,
  type InventoryStockTableSorting,
} from '~/types/inventory'
import type { PaginatedResponse } from '~/types/api'

export type InventoryStockPageQueryParams = {
  preset: InventoryStockPreset
  page: number
  pageSize: number
  search: string
  sorting: InventoryStockTableSorting[]
  deviceTypeId?: number | null
}

const getStockEndpointForPreset = (preset: InventoryStockPreset): string => {
  if (preset === 'favorite') {
    return `${INVENTORY_STOCKS_ENDPOINT}favorites/`
  }

  if (preset === 'low_stock') {
    return `${INVENTORY_STOCKS_ENDPOINT}low_stock/`
  }

  if (preset === 'expired') {
    return `${INVENTORY_STOCKS_ENDPOINT}expired/`
  }

  if (preset === 'archived') {
    return `${INVENTORY_STOCKS_ENDPOINT}archived/`
  }

  return INVENTORY_STOCKS_ENDPOINT
}

const toOrderingValue = (sorting: InventoryStockTableSorting[]): string => {
  return sorting.map((sortRule) => (sortRule.desc ? `-${sortRule.id}` : sortRule.id)).join(',')
}

const fetchInventoryStockPage = async (
  params: InventoryStockPageQueryParams,
): Promise<PaginatedResponse<InventoryStockListItem>> => {
  const endpoint = getStockEndpointForPreset(params.preset)
  const requestParams: Record<string, string> = {
    page: String(params.page),
    page_size: String(params.pageSize),
  }

  if (params.search.trim() !== '') {
    requestParams.search = params.search.trim()
  }

  if (params.deviceTypeId && params.deviceTypeId > 0) {
    requestParams.device_type = String(params.deviceTypeId)
  }

  const orderingValue = toOrderingValue(params.sorting)
  if (orderingValue !== '') {
    requestParams.ordering = orderingValue
  }

  const { data, error } = await useAPI<PaginatedResponse<InventoryStockListItem>>(endpoint, {
    method: 'GET',
    params: requestParams,
  })

  if (error.value || !data.value) {
    throw (error.value ?? new Error(INVENTORY_STOCKS_ERROR_MESSAGE)) as Error
  }

  return data.value
}

export const useInventoryStockPageQuery = (
  paramsRef: ComputedRef<InventoryStockPageQueryParams>,
  enabledRef?: ComputedRef<boolean>,
) =>
  useQuery({
    queryKey: computed(() => ['inventory-stocks-page', paramsRef.value]),
    enabled: computed(() => enabledRef?.value ?? true),
    queryFn: () => fetchInventoryStockPage(paramsRef.value),
  })
