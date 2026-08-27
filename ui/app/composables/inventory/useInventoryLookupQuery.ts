import { computed, type ComputedRef } from 'vue'
import { useQuery } from '@tanstack/vue-query'
import { useAPI } from '~/composables/useAPI'
import { parseListPage, type ListResponse } from '~/composables/inventory/listResponse'
import {
  INVENTORY_ATTRIBUTES_ENDPOINT,
  INVENTORY_BRANDS_ENDPOINT,
  INVENTORY_DEVICE_TYPES_ENDPOINT,
  INVENTORY_ITEM_TYPES_ENDPOINT,
  INVENTORY_LOOKUPS_QUERY_KEY,
  INVENTORY_MANUFACTURERS_ENDPOINT,
  INVENTORY_ROOMS_ENDPOINT,
  INVENTORY_SECTORS_ENDPOINT,
  INVENTORY_UNITS_ENDPOINT,
  INVENTORY_VENDORS_ENDPOINT,
  type InventoryLookupItem,
  type InventoryRoom,
  type InventorySector,
} from '~/types/inventory'

const fetchAllPages = async <T>(endpoint: string): Promise<T[]> => {
  const allItems: T[] = []
  const visitedEndpoints = new Set<string>()
  let nextEndpoint: string | null = endpoint

  while (nextEndpoint) {
    if (visitedEndpoints.has(nextEndpoint)) break
    visitedEndpoints.add(nextEndpoint)

    const { data, error } = await useAPI<ListResponse<T>>(nextEndpoint, { method: 'GET' })

    if (error.value || !data.value) {
      throw (error.value ?? new Error(`Failed to load ${endpoint}`)) as Error
    }

    const parsedPage = parseListPage(data.value)
    allItems.push(...parsedPage.items)
    nextEndpoint = parsedPage.nextEndpoint
  }

  return allItems
}

export type InventoryLookupBundle = {
  rooms: InventoryRoom[]
  sectors: InventorySector[]
  manufacturers: InventoryLookupItem[]
  brands: InventoryLookupItem[]
  vendors: InventoryLookupItem[]
  deviceTypes: InventoryLookupItem[]
  itemTypes: InventoryLookupItem[]
  attributes: InventoryLookupItem[]
  units: InventoryLookupItem[]
}

const fetchInventoryLookups = async (): Promise<InventoryLookupBundle> => {
  const [rooms, sectors, manufacturers, brands, vendors, deviceTypes, itemTypes, attributes, units] = await Promise.all(
    [
      fetchAllPages<InventoryRoom>(INVENTORY_ROOMS_ENDPOINT),
      fetchAllPages<InventorySector>(INVENTORY_SECTORS_ENDPOINT),
      fetchAllPages<InventoryLookupItem>(INVENTORY_MANUFACTURERS_ENDPOINT),
      fetchAllPages<InventoryLookupItem>(INVENTORY_BRANDS_ENDPOINT),
      fetchAllPages<InventoryLookupItem>(INVENTORY_VENDORS_ENDPOINT),
      fetchAllPages<InventoryLookupItem>(INVENTORY_DEVICE_TYPES_ENDPOINT),
      fetchAllPages<InventoryLookupItem>(INVENTORY_ITEM_TYPES_ENDPOINT),
      fetchAllPages<InventoryLookupItem>(INVENTORY_ATTRIBUTES_ENDPOINT),
      fetchAllPages<InventoryLookupItem>(INVENTORY_UNITS_ENDPOINT),
    ],
  )

  return {
    rooms,
    sectors,
    manufacturers,
    brands,
    vendors,
    deviceTypes,
    itemTypes,
    attributes,
    units,
  }
}

export const useInventoryLookupsQuery = (enabledRef?: ComputedRef<boolean>) =>
  useQuery({
    queryKey: INVENTORY_LOOKUPS_QUERY_KEY,
    enabled: computed(() => enabledRef?.value ?? true),
    queryFn: fetchInventoryLookups,
  })
