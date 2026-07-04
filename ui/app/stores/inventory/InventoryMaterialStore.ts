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

  const buildMaterialRequestBody = (
    payload: CreateInventoryMaterialPayload | UpdateInventoryMaterialPayload,
  ): CreateInventoryMaterialPayload | UpdateInventoryMaterialPayload | FormData => {
    const hasSafetyDataSheetFile =
      typeof File !== 'undefined' && payload.safety_data_sheet instanceof File

    if (!hasSafetyDataSheetFile) {
      return payload
    }

    const formData = new FormData()

    Object.entries(payload).forEach(([key, value]) => {
      if (value == null) {
        return
      }

      if (Array.isArray(value)) {
        value.forEach((entryValue) => {
          formData.append(key, String(entryValue))
        })
        return
      }

      if (typeof File !== 'undefined' && value instanceof File) {
        formData.append(key, value)
        return
      }

      formData.append(key, String(value))
    })

    return formData
  }

  const createMaterial = async (payload: CreateInventoryMaterialPayload): Promise<InventoryMaterialDetail> => {
    isCreatingMaterial.value = true
    try {
      return await requestApiData<InventoryMaterialDetail>(
        'inventory/materials/',
        {
          method: 'POST',
          body: buildMaterialRequestBody(payload),
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
          body: buildMaterialRequestBody(payload),
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
