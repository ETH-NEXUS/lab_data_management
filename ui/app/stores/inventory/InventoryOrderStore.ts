import { defineStore } from 'pinia'
import { ref } from 'vue'
import {
  INVENTORY_CREATE_ORDER_ERROR_MESSAGE,
  INVENTORY_UPDATE_ORDER_ERROR_MESSAGE,
  type CreateInventoryOrderPayload,
  type InventoryOrderDetail,
  type UpdateInventoryOrderPayload,
} from '~/types/inventory'
import { requestApiData } from '~/utils/apiRequests'

export const useInventoryOrderStore = defineStore('inventoryOrderStore', () => {
  const isCreatingOrder = ref(false)
  const isUpdatingOrder = ref(false)

  const createOrder = async (payload: CreateInventoryOrderPayload): Promise<InventoryOrderDetail> => {
    isCreatingOrder.value = true
    try {
      return await requestApiData<InventoryOrderDetail>(
        'inventory/orders/',
        {
          method: 'POST',
          body: payload,
        },
        INVENTORY_CREATE_ORDER_ERROR_MESSAGE,
      )
    } finally {
      isCreatingOrder.value = false
    }
  }

  const updateOrder = async (orderId: number, payload: UpdateInventoryOrderPayload): Promise<InventoryOrderDetail> => {
    isUpdatingOrder.value = true
    try {
      return await requestApiData<InventoryOrderDetail>(
        `inventory/orders/${orderId}/`,
        {
          method: 'PATCH',
          body: payload,
        },
        INVENTORY_UPDATE_ORDER_ERROR_MESSAGE,
      )
    } finally {
      isUpdatingOrder.value = false
    }
  }

  return {
    isCreatingOrder,
    isUpdatingOrder,
    createOrder,
    updateOrder,
  }
})
