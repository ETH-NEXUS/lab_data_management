import { defineStore } from 'pinia'
import { ref } from 'vue'
import {
  INVENTORY_ARCHIVE_STOCK_ERROR_MESSAGE,
  INVENTORY_CREATE_STOCK_ERROR_MESSAGE,
  INVENTORY_MARK_FAVORITE_ERROR_MESSAGE,
  INVENTORY_RESTORE_STOCK_ERROR_MESSAGE,
  INVENTORY_UNMARK_FAVORITE_ERROR_MESSAGE,
  INVENTORY_UPDATE_STOCK_ERROR_MESSAGE,
  type CreateInventoryStockPayload,
  type InventoryStockDetail,
  type UpdateInventoryStockPayload,
} from '~/types/inventory'
import { requestApiData } from '~/utils/apiRequests'

export const useInventoryStockStore = defineStore('inventoryStockStore', () => {
  const isCreatingStock = ref(false)
  const isUpdatingStock = ref(false)
  const isMarkingFavorite = ref(false)
  const isUnmarkingFavorite = ref(false)
  const isArchivingStock = ref(false)
  const isRestoringStock = ref(false)

  const createStock = async (payload: CreateInventoryStockPayload): Promise<InventoryStockDetail> => {
    isCreatingStock.value = true
    try {
      return await requestApiData<InventoryStockDetail>(
        'inventory/stocks/',
        {
          method: 'POST',
          body: payload,
        },
        INVENTORY_CREATE_STOCK_ERROR_MESSAGE,
      )
    } finally {
      isCreatingStock.value = false
    }
  }

  const updateStock = async (stockId: number, payload: UpdateInventoryStockPayload): Promise<InventoryStockDetail> => {
    isUpdatingStock.value = true
    try {
      return await requestApiData<InventoryStockDetail>(
        `inventory/stocks/${stockId}/`,
        {
          method: 'PATCH',
          body: payload,
        },
        INVENTORY_UPDATE_STOCK_ERROR_MESSAGE,
      )
    } finally {
      isUpdatingStock.value = false
    }
  }

  const markFavorite = async (stockId: number): Promise<InventoryStockDetail> => {
    isMarkingFavorite.value = true
    try {
      return await requestApiData<InventoryStockDetail>(
        `inventory/stocks/${stockId}/mark_favorite/`,
        {
          method: 'POST',
        },
        INVENTORY_MARK_FAVORITE_ERROR_MESSAGE,
      )
    } finally {
      isMarkingFavorite.value = false
    }
  }

  const unmarkFavorite = async (stockId: number): Promise<InventoryStockDetail> => {
    isUnmarkingFavorite.value = true
    try {
      return await requestApiData<InventoryStockDetail>(
        `inventory/stocks/${stockId}/unmark_favorite/`,
        {
          method: 'POST',
        },
        INVENTORY_UNMARK_FAVORITE_ERROR_MESSAGE,
      )
    } finally {
      isUnmarkingFavorite.value = false
    }
  }

  const archiveStock = async (stockId: number): Promise<InventoryStockDetail> => {
    isArchivingStock.value = true
    try {
      return await requestApiData<InventoryStockDetail>(
        `inventory/stocks/${stockId}/archive/`,
        {
          method: 'POST',
        },
        INVENTORY_ARCHIVE_STOCK_ERROR_MESSAGE,
      )
    } finally {
      isArchivingStock.value = false
    }
  }

  const restoreStock = async (stockId: number): Promise<InventoryStockDetail> => {
    isRestoringStock.value = true
    try {
      return await requestApiData<InventoryStockDetail>(
        `inventory/stocks/${stockId}/restore/`,
        {
          method: 'POST',
        },
        INVENTORY_RESTORE_STOCK_ERROR_MESSAGE,
      )
    } finally {
      isRestoringStock.value = false
    }
  }

  return {
    isCreatingStock,
    isUpdatingStock,
    isMarkingFavorite,
    isUnmarkingFavorite,
    isArchivingStock,
    isRestoringStock,
    createStock,
    updateStock,
    markFavorite,
    unmarkFavorite,
    archiveStock,
    restoreStock,
  }
})
