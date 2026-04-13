import { computed, type ComputedRef } from 'vue'
import { useQuery } from '@tanstack/vue-query'
import { useAPI } from '~/composables/useAPI'
import { parseListPage, type ListResponse } from '~/composables/inventory/listResponse'
import {
  INVENTORY_MATERIALS_ENDPOINT,
  INVENTORY_MATERIALS_ERROR_MESSAGE,
  INVENTORY_MATERIAL_ERROR_MESSAGE,
  INVENTORY_MATERIALS_QUERY_KEY,
  getInventoryMaterialQueryKey,
  type InventoryMaterialDetail,
  type InventoryMaterialListItem,
} from '~/types/inventory'

const fetchInventoryMaterials = async (): Promise<InventoryMaterialListItem[]> => {
  const allItems: InventoryMaterialListItem[] = []
  const visitedEndpoints = new Set<string>()
  let nextEndpoint: string | null = INVENTORY_MATERIALS_ENDPOINT

  while (nextEndpoint) {
    if (visitedEndpoints.has(nextEndpoint)) break
    visitedEndpoints.add(nextEndpoint)

    const { data, error } = await useAPI<ListResponse<InventoryMaterialListItem>>(nextEndpoint, {
      method: 'GET',
    })

    if (error.value || !data.value) {
      throw (error.value ?? new Error(INVENTORY_MATERIALS_ERROR_MESSAGE)) as Error
    }

    const parsedPage = parseListPage(data.value)
    allItems.push(...parsedPage.items)
    nextEndpoint = parsedPage.nextEndpoint
  }

  return allItems
}

const fetchInventoryMaterial = async (materialId: number): Promise<InventoryMaterialDetail> => {
  const { data, error } = await useAPI<InventoryMaterialDetail>(`inventory/materials/${materialId}/`, {
    method: 'GET',
  })

  if (error.value || !data.value) {
    throw (error.value ?? new Error(INVENTORY_MATERIAL_ERROR_MESSAGE)) as Error
  }

  return data.value
}

export const useInventoryMaterialsQuery = () =>
  useQuery({
    queryKey: INVENTORY_MATERIALS_QUERY_KEY,
    queryFn: fetchInventoryMaterials,
  })

export const useInventoryMaterialQuery = (materialIdRef: ComputedRef<number>) =>
  useQuery({
    queryKey: computed(() => getInventoryMaterialQueryKey(materialIdRef.value)),
    enabled: computed(() => Number.isInteger(materialIdRef.value) && materialIdRef.value > 0),
    queryFn: () => fetchInventoryMaterial(materialIdRef.value),
  })
