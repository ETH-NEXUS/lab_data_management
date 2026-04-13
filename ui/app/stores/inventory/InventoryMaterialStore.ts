import { defineStore } from 'pinia'
import { ref } from 'vue'
import {
  INVENTORY_CREATE_MATERIAL_ERROR_MESSAGE,
  INVENTORY_UPDATE_MATERIAL_ERROR_MESSAGE,
  type CreateInventoryMaterialPayload,
  type InventoryMaterialDetail,
  type UpdateInventoryMaterialPayload,
} from '~/types/inventory'
import { requestApiData } from '~/utils/apiRequests'

export const useInventoryMaterialStore = defineStore('inventoryMaterialStore', () => {
  const isCreatingMaterial = ref(false)
  const isUpdatingMaterial = ref(false)

  const createMaterial = async (payload: CreateInventoryMaterialPayload): Promise<InventoryMaterialDetail> => {
    isCreatingMaterial.value = true
    try {
      return await requestApiData<InventoryMaterialDetail>(
        'inventory/materials/',
        {
          method: 'POST',
          body: payload,
        },
        INVENTORY_CREATE_MATERIAL_ERROR_MESSAGE,
      )
    } finally {
      isCreatingMaterial.value = false
    }
  }

  const updateMaterial = async (
    materialId: number,
    payload: UpdateInventoryMaterialPayload,
  ): Promise<InventoryMaterialDetail> => {
    isUpdatingMaterial.value = true
    try {
      return await requestApiData<InventoryMaterialDetail>(
        `inventory/materials/${materialId}/`,
        {
          method: 'PATCH',
          body: payload,
        },
        INVENTORY_UPDATE_MATERIAL_ERROR_MESSAGE,
      )
    } finally {
      isUpdatingMaterial.value = false
    }
  }

  return {
    isCreatingMaterial,
    isUpdatingMaterial,
    createMaterial,
    updateMaterial,
  }
})
