import { defineStore } from 'pinia'
import { computed, ref, watch } from 'vue'
import type { ColumnFiltersState, ColumnOrderState, SortingState, VisibilityState } from '@tanstack/vue-table'
import { useInventoryStockTablePreferenceQuery } from '~/composables/inventory/useInventoryStockTablePreferenceQuery'
import { requestApiData } from '~/utils/apiRequests'
import {
  INVENTORY_STOCK_TABLE_PREFERENCES_ENDPOINT,
  INVENTORY_UPDATE_STOCK_TABLE_PREFERENCE_ERROR_MESSAGE,
  type InventoryStockPreset,
  type InventoryStockTablePreference,
} from '~/types/inventory'

/**
 * Provides normalized inventory stock table preference state from the backend.
 *
 * This store intentionally starts read-only. Save-back behavior will be added
 * in a later step once the consuming UI is fully rewired.
 */
export const useInventoryStockTablePreferenceStore = defineStore('inventoryStockTablePreferenceStore', () => {
  const preferenceQuery = useInventoryStockTablePreferenceQuery()
  const isSavingPreference = ref(false)

  const presetState = ref<InventoryStockPreset>('all')
  const sortingState = ref<SortingState>([])
  const globalFilterState = ref('')
  const columnFiltersState = ref<ColumnFiltersState>([])
  const columnOrderState = ref<ColumnOrderState>([])
  const columnVisibilityState = ref<VisibilityState>({})

  const preference = computed<InventoryStockTablePreference | null>(() => {
    return preferenceQuery.data.value ?? null
  })

  watch(
    () => preference.value,
    (nextPreference) => {
      presetState.value = nextPreference?.preset ?? 'all'
      sortingState.value = nextPreference?.sorting ?? []
      globalFilterState.value = ''
      columnFiltersState.value = nextPreference?.column_filters ?? []
      columnOrderState.value = nextPreference?.column_order ?? []
      columnVisibilityState.value = nextPreference?.column_visibility ?? {}
    },
    { immediate: true },
  )

  const savePreference = async (payload: Partial<InventoryStockTablePreference>): Promise<void> => {
    isSavingPreference.value = true
    try {
      await requestApiData<InventoryStockTablePreference>(
        `${INVENTORY_STOCK_TABLE_PREFERENCES_ENDPOINT}current/`,
        {
          method: 'PATCH',
          body: payload,
        },
        INVENTORY_UPDATE_STOCK_TABLE_PREFERENCE_ERROR_MESSAGE,
      )
    } finally {
      isSavingPreference.value = false
    }
  }

  const updateSortingState = async (nextSortingState: SortingState): Promise<void> => {
    sortingState.value = nextSortingState
    await savePreference({ sorting: nextSortingState })
  }

  const updatePresetState = async (nextPresetState: InventoryStockPreset): Promise<void> => {
    presetState.value = nextPresetState
    await savePreference({ preset: nextPresetState })
  }

  const updateGlobalFilterState = async (nextGlobalFilterState: string): Promise<void> => {
    globalFilterState.value = nextGlobalFilterState
  }

  const updateColumnFiltersState = async (nextColumnFiltersState: ColumnFiltersState): Promise<void> => {
    columnFiltersState.value = nextColumnFiltersState
    await savePreference({ column_filters: nextColumnFiltersState })
  }

  const updateColumnOrderState = async (nextColumnOrderState: ColumnOrderState): Promise<void> => {
    columnOrderState.value = nextColumnOrderState
    await savePreference({ column_order: nextColumnOrderState })
  }

  const updateColumnVisibilityState = async (nextColumnVisibilityState: VisibilityState): Promise<void> => {
    columnVisibilityState.value = nextColumnVisibilityState
    await savePreference({ column_visibility: nextColumnVisibilityState })
  }

  return {
    preferenceQuery,
    preference,
    isSavingPreference,
    presetState,
    sortingState,
    globalFilterState,
    columnFiltersState,
    columnOrderState,
    columnVisibilityState,
    updatePresetState,
    updateSortingState,
    updateGlobalFilterState,
    updateColumnFiltersState,
    updateColumnOrderState,
    updateColumnVisibilityState,
  }
})
