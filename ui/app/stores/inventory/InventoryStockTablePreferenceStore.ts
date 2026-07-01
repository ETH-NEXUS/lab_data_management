import { defineStore } from 'pinia'
import { computed } from 'vue'
import type { ColumnFiltersState, ColumnOrderState, SortingState, VisibilityState } from '@tanstack/vue-table'
import { useInventoryStockTablePreferenceQuery } from '~/composables/inventory/useInventoryStockTablePreferenceQuery'
import type { InventoryStockTablePreference } from '~/types/inventory'

/**
 * Provides normalized inventory stock table preference state from the backend.
 *
 * This store intentionally starts read-only. Save-back behavior will be added
 * in a later step once the consuming UI is fully rewired.
 */
export const useInventoryStockTablePreferenceStore = defineStore('inventoryStockTablePreferenceStore', () => {
  const preferenceQuery = useInventoryStockTablePreferenceQuery()

  const preference = computed<InventoryStockTablePreference | null>(() => {
    return preferenceQuery.data.value ?? null
  })

  const sortingState = computed<SortingState>(() => {
    return preference.value?.sorting ?? []
  })

  const globalFilterState = computed<string>(() => {
    return ''
  })

  const columnFiltersState = computed<ColumnFiltersState>(() => {
    return preference.value?.column_filters ?? []
  })

  const columnOrderState = computed<ColumnOrderState>(() => {
    return preference.value?.column_order ?? []
  })

  const columnVisibilityState = computed<VisibilityState>(() => {
    return preference.value?.column_visibility ?? {}
  })

  return {
    preferenceQuery,
    preference,
    sortingState,
    globalFilterState,
    columnFiltersState,
    columnOrderState,
    columnVisibilityState,
  }
})
