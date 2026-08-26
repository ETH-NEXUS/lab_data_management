import { computed, type ComputedRef } from 'vue'
import { useQuery } from '@tanstack/vue-query'
import { useAPI } from '~/composables/useAPI'
import { parseListPage, type ListResponse } from '~/composables/inventory/listResponse'
import {
  INVENTORY_ORDERS_ENDPOINT,
  INVENTORY_ORDERS_ERROR_MESSAGE,
  INVENTORY_ORDER_ERROR_MESSAGE,
  INVENTORY_ORDERS_QUERY_KEY,
  getInventoryOrderQueryKey,
  type InventoryOrderDetail,
  type InventoryOrderListItem,
} from '~/types/inventory'

const fetchInventoryOrders = async (): Promise<InventoryOrderListItem[]> => {
  const allItems: InventoryOrderListItem[] = []
  const visitedEndpoints = new Set<string>()
  let nextEndpoint: string | null = INVENTORY_ORDERS_ENDPOINT

  while (nextEndpoint) {
    if (visitedEndpoints.has(nextEndpoint)) break
    visitedEndpoints.add(nextEndpoint)

    const { data, error } = await useAPI<ListResponse<InventoryOrderListItem>>(nextEndpoint, {
      method: 'GET',
    })

    if (error.value || !data.value) {
      throw (error.value ?? new Error(INVENTORY_ORDERS_ERROR_MESSAGE)) as Error
    }

    const parsedPage = parseListPage(data.value)
    allItems.push(...parsedPage.items)
    nextEndpoint = parsedPage.nextEndpoint
  }

  return allItems
}

const fetchInventoryOrder = async (orderId: number): Promise<InventoryOrderDetail> => {
  const { data, error } = await useAPI<InventoryOrderDetail>(`inventory/orders/${orderId}/`, {
    method: 'GET',
  })

  if (error.value || !data.value) {
    throw (error.value ?? new Error(INVENTORY_ORDER_ERROR_MESSAGE)) as Error
  }

  return data.value
}

const fetchAwaitingCheckInOrders = async (): Promise<InventoryOrderListItem[]> => {
  const { data, error } = await useAPI<InventoryOrderListItem[]>(`${INVENTORY_ORDERS_ENDPOINT}awaiting_check_in/`, {
    method: 'GET',
  })

  if (error.value || !data.value) {
    throw (error.value ?? new Error(INVENTORY_ORDERS_ERROR_MESSAGE)) as Error
  }

  return data.value
}

export const useInventoryOrdersQuery = () =>
  useQuery({
    queryKey: INVENTORY_ORDERS_QUERY_KEY,
    queryFn: fetchInventoryOrders,
  })

export const useInventoryAwaitingCheckInOrdersQuery = () =>
  useQuery({
    queryKey: [...INVENTORY_ORDERS_QUERY_KEY, 'awaiting-check-in'],
    queryFn: fetchAwaitingCheckInOrders,
  })

export const useInventoryOrderQuery = (orderIdRef: ComputedRef<number>) =>
  useQuery({
    queryKey: computed(() => getInventoryOrderQueryKey(orderIdRef.value)),
    enabled: computed(() => Number.isInteger(orderIdRef.value) && orderIdRef.value > 0),
    queryFn: () => fetchInventoryOrder(orderIdRef.value),
  })
