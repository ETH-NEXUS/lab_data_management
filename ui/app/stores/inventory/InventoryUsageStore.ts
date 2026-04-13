import { defineStore } from 'pinia'
import { ref } from 'vue'
import {
  INVENTORY_CREATE_USAGE_ERROR_MESSAGE,
  INVENTORY_UPDATE_USAGE_ERROR_MESSAGE,
  type CreateInventoryUsagePayload,
  type InventoryUsageDetail,
  type UpdateInventoryUsagePayload,
} from '~/types/inventory'
import { requestApiData } from '~/utils/apiRequests'

export const useInventoryUsageStore = defineStore('inventoryUsageStore', () => {
  const isCreatingUsage = ref(false)
  const isUpdatingUsage = ref(false)

  const createUsage = async (payload: CreateInventoryUsagePayload): Promise<InventoryUsageDetail> => {
    isCreatingUsage.value = true
    try {
      return await requestApiData<InventoryUsageDetail>(
        'inventory/material-usages/',
        {
          method: 'POST',
          body: payload,
        },
        INVENTORY_CREATE_USAGE_ERROR_MESSAGE,
      )
    } finally {
      isCreatingUsage.value = false
    }
  }

  const updateUsage = async (usageId: number, payload: UpdateInventoryUsagePayload): Promise<InventoryUsageDetail> => {
    isUpdatingUsage.value = true
    try {
      return await requestApiData<InventoryUsageDetail>(
        `inventory/material-usages/${usageId}/`,
        {
          method: 'PATCH',
          body: payload,
        },
        INVENTORY_UPDATE_USAGE_ERROR_MESSAGE,
      )
    } finally {
      isUpdatingUsage.value = false
    }
  }

  return {
    isCreatingUsage,
    isUpdatingUsage,
    createUsage,
    updateUsage,
  }
})
