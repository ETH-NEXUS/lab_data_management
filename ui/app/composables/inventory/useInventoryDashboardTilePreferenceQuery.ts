import { ref } from 'vue'
import { useQuery, useQueryClient } from '@tanstack/vue-query'
import { useAPI } from '~/composables/useAPI'
import {
  INVENTORY_DASHBOARD_TILE_PREFERENCES_ENDPOINT,
  INVENTORY_DASHBOARD_TILE_PREFERENCES_QUERY_KEY,
  INVENTORY_DASHBOARD_TILE_PREFERENCE_ERROR_MESSAGE,
  INVENTORY_UPDATE_DASHBOARD_TILE_PREFERENCE_ERROR_MESSAGE,
  type InventoryDashboardTilePreference,
  type UpdateInventoryDashboardTilePreferencePayload,
} from '~/types/inventory'

const CURRENT_DASHBOARD_TILE_PREFERENCES_ENDPOINT = `${INVENTORY_DASHBOARD_TILE_PREFERENCES_ENDPOINT}current/`

const fetchDashboardTilePreferences = async (): Promise<InventoryDashboardTilePreference[]> => {
  const { data, error } = await useAPI<InventoryDashboardTilePreference[]>(
    CURRENT_DASHBOARD_TILE_PREFERENCES_ENDPOINT,
    {
      method: 'GET',
    },
  )

  if (error.value || !data.value) {
    throw (error.value ?? new Error(INVENTORY_DASHBOARD_TILE_PREFERENCE_ERROR_MESSAGE)) as Error
  }

  return data.value
}

const updateDashboardTilePreferences = async (
  payload: UpdateInventoryDashboardTilePreferencePayload,
): Promise<InventoryDashboardTilePreference[]> => {
  const { data, error } = await useAPI<InventoryDashboardTilePreference[]>(
    CURRENT_DASHBOARD_TILE_PREFERENCES_ENDPOINT,
    {
      method: 'PUT',
      body: payload,
    },
  )

  if (error.value || !data.value) {
    throw (error.value ?? new Error(INVENTORY_UPDATE_DASHBOARD_TILE_PREFERENCE_ERROR_MESSAGE)) as Error
  }

  return data.value
}

export const useInventoryDashboardTilePreferenceQuery = () =>
  useQuery({
    queryKey: INVENTORY_DASHBOARD_TILE_PREFERENCES_QUERY_KEY,
    queryFn: fetchDashboardTilePreferences,
  })

export const useInventoryDashboardTilePreferences = () => {
  const queryClient = useQueryClient()
  const preferenceQuery = useInventoryDashboardTilePreferenceQuery()
  const isSaving = ref(false)

  const saveTileKeys = async (tileKeys: string[]): Promise<void> => {
    isSaving.value = true

    try {
      const preferences = await updateDashboardTilePreferences({ tile_keys: tileKeys })
      queryClient.setQueryData(INVENTORY_DASHBOARD_TILE_PREFERENCES_QUERY_KEY, preferences)
    } finally {
      isSaving.value = false
    }
  }

  return {
    preferenceQuery,
    isSaving,
    saveTileKeys,
  }
}
