import { computed } from 'vue'
import { useQuery } from '@tanstack/vue-query'
import { useAPI } from '~/composables/useAPI'
import {
  INVENTORY_STOCK_TABLE_PREFERENCES_ENDPOINT,
  INVENTORY_STOCK_TABLE_PREFERENCE_ERROR_MESSAGE,
  getInventoryStockTablePreferenceQueryKey,
  type InventoryStockTablePreference,
} from '~/types/inventory'

const TABLE_KEY_INVENTORY_STOCK = 'inventory_stock'

const fetchInventoryStockTablePreference = async (): Promise<InventoryStockTablePreference> => {
  const { data, error } = await useAPI<InventoryStockTablePreference>(
    `${INVENTORY_STOCK_TABLE_PREFERENCES_ENDPOINT}current/`,
    {
      method: 'GET',
    },
  )

  if (error.value || !data.value) {
    throw (error.value ?? new Error(INVENTORY_STOCK_TABLE_PREFERENCE_ERROR_MESSAGE)) as Error
  }

  return data.value
}

export const useInventoryStockTablePreferenceQuery = () =>
  useQuery({
    queryKey: computed(() => getInventoryStockTablePreferenceQueryKey(TABLE_KEY_INVENTORY_STOCK)),
    queryFn: fetchInventoryStockTablePreference,
  })
