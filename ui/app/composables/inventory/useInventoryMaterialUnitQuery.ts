import { computed, type ComputedRef } from 'vue'
import { useQuery } from '@tanstack/vue-query'
import { useAPI } from '~/composables/useAPI'
import { parseListPage, type ListResponse } from '~/composables/inventory/listResponse'
import type { InventoryMaterialUnitSummary } from '~/types/inventory'

const INVENTORY_MATERIAL_UNITS_ERROR_MESSAGE = 'Failed to load inventory material units.'
const INVENTORY_MATERIAL_UNITS_QUERY_KEY = ['inventory-material-units']

const fetchInventoryMaterialOrderUnits = async (materialId: number): Promise<InventoryMaterialUnitSummary[]> => {
  const allItems: InventoryMaterialUnitSummary[] = []
  const visitedEndpoints = new Set<string>()
  let nextEndpoint: string | null = `inventory/material-units/?material=${materialId}&is_order_unit=true`

  while (nextEndpoint) {
    if (visitedEndpoints.has(nextEndpoint)) break
    visitedEndpoints.add(nextEndpoint)

    const { data, error } = await useAPI<ListResponse<InventoryMaterialUnitSummary>>(nextEndpoint, {
      method: 'GET',
    })

    if (error.value || !data.value) {
      throw (error.value ?? new Error(INVENTORY_MATERIAL_UNITS_ERROR_MESSAGE)) as Error
    }

    const parsedPage = parseListPage(data.value)
    allItems.push(...parsedPage.items)
    nextEndpoint = parsedPage.nextEndpoint
  }

  return allItems
}

export const useInventoryMaterialOrderUnitsQuery = (materialIdRef: ComputedRef<number>) =>
  useQuery({
    queryKey: computed(() => [...INVENTORY_MATERIAL_UNITS_QUERY_KEY, materialIdRef.value, 'order']),
    enabled: computed(() => Number.isInteger(materialIdRef.value) && materialIdRef.value > 0),
    queryFn: () => fetchInventoryMaterialOrderUnits(materialIdRef.value),
  })
