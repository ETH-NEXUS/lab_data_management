import { computed, nextTick, type ComputedRef } from 'vue'
import { useQueryClient } from '@tanstack/vue-query'
import { useInventoryMaterialStore } from '~/stores/inventory/InventoryMaterialStore'
import { useInventoryStockStore } from '~/stores/inventory/InventoryStock'
import {
  getInventoryMaterialQueryKey,
  INVENTORY_MATERIALS_QUERY_KEY,
  INVENTORY_STOCK_PAGES_QUERY_KEY,
  INVENTORY_STOCKS_QUERY_KEY,
  type CreateInventoryStockPayload,
  type UpdateInventoryMaterialPayload,
} from '~/types/inventory'
import {
  formatDecimal,
  parseDecimal,
  parseInteger,
} from '~/components/inventory/add-item/inventoryAddItemForm.utils'
import { parsePositiveIntegerList } from '~/components/inventory/inventorySectorSelection.utils'
import { useInventoryAddItemFormState } from '~/components/inventory/add-item/useInventoryAddItemFormState'
import { getErrorMessage } from '~/utils/errors'

type UseInventoryAddItemFormParams = {
  open: ComputedRef<boolean>
  onSaved: () => void
}

/**
 * Manages add-item modal form state, dependent lookup data, validation, and stock creation submit.
 */
export const useInventoryAddItemForm = ({ open, onSaved }: UseInventoryAddItemFormParams) => {
  const toast = useToast()
  const queryClient = useQueryClient()
  const inventoryMaterialStore = useInventoryMaterialStore()
  const inventoryStockStore = useInventoryStockStore()
  const isSubmitting = computed<boolean>(() => {
    return inventoryStockStore.isCreatingStock || inventoryMaterialStore.isUpdatingMaterial
  })

  const {
    formState,
    selectedOrderId,
    selectedItemTypeId,
    selectedMaterialId,
    selectedRoomId,
    roomOptions,
    selectedMaterialQuery,
    sortedMaterials,
    filteredMaterials,
    sortedOrders,
    filteredOrders,
    selectedOrder,
    orderPrefillOptions,
    sectorOptions,
    stockUnitOptions,
    isStockUnitsLoading,
    materialsQuery,
    lookupsQuery,
    ordersQuery,
    materialsErrorMessage,
    sectorsErrorMessage,
    stockUnitsErrorMessage,
    ordersErrorMessage,
  } = useInventoryAddItemFormState({ open })

  /**
   * Builds read-only supplier and catalog display values from the selected material.
   *
   * Returned data example:
   * - `{ manufacturer: 'Corning', vendor: 'Huberlab', manufacturerCatalogNumber: '3005', vendorCatalogNumber: 'H-42' }`
   * - `{ manufacturer: null, vendor: null, manufacturerCatalogNumber: null, vendorCatalogNumber: null }`
   */
  const supplierDetails = computed<{
    manufacturer: string | null
    vendor: string | null
    manufacturerCatalogNumber: string | null
    vendorCatalogNumber: string | null
  }>(() => {
    const material = selectedMaterialQuery.data.value

    return {
      manufacturer: material?.manufacturer?.label || material?.manufacturer?.name || null,
      vendor: material?.vendor?.label || material?.vendor?.name || null,
      manufacturerCatalogNumber: material?.manufacturer_catalog_number || null,
      vendorCatalogNumber: material?.vendor_catalog_number || null,
    }
  })

  const isSelectedMaterialReagent = computed<boolean>(() => {
    const itemTypeName =
      selectedMaterialQuery.data.value?.item_type?.name ||
      selectedMaterialQuery.data.value?.item_type?.label ||
      ''

    return itemTypeName.trim().toLowerCase() === 'reagent'
  })

  const validationMessages = computed<string[]>(() => {
    const messages: string[] = []
    const itemTypeId = Number.parseInt(formState.value.itemTypeId, 10)
    const materialId = Number.parseInt(formState.value.materialId, 10)
    const sectorIds = parsePositiveIntegerList(formState.value.sectorIds)
    const stockUnitId = Number.parseInt(formState.value.stockUnitId, 10)
    const quantityValue = parseDecimal(formState.value.quantity)
    const minimumQuantityValue = parseInteger(formState.value.minimumQuantity)

    if (!Number.isInteger(itemTypeId) || itemTypeId <= 0) {
      messages.push('Item type is required before selecting a material or order.')
    }
    if (!Number.isInteger(materialId) || materialId <= 0) {
      messages.push('Material is required.')
    }
    if (isSelectedMaterialReagent.value && formState.value.reagentStorageTemperature.trim() === '') {
      messages.push('Storage temperature is required for reagents.')
    }
    if (sectorIds.length === 0) {
      messages.push('At least one sector is required.')
    }
    if (!Number.isInteger(stockUnitId) || stockUnitId <= 0) {
      messages.push('Stock unit is required.')
    }
    if (quantityValue === null || quantityValue <= 0) {
      messages.push('Quantity must be greater than zero.')
    }
    if (minimumQuantityValue === null || minimumQuantityValue < 0) {
      messages.push('Minimum quantity must be a whole number that is zero or greater.')
    }

    return messages
  })

  const isFormValid = computed<boolean>(() => validationMessages.value.length === 0)

  const isOrderPrefillDisabled = computed<boolean>(() => {
    if (isSubmitting.value) {
      return true
    }

    const orderId = Number.parseInt(selectedOrderId.value, 10)
    return !Number.isInteger(orderId) || orderId <= 0
  })

  const isSaveDisabled = computed<boolean>(() => {
    if (isSubmitting.value) {
      return true
    }
    return !isFormValid.value
  })

  /**
   * Builds API payload for inventory stock creation from form draft state.
   *
   * Accepted draft example:
   * - `{ materialId: '12', sectorIds: ['7', '8'], stockUnitId: '31', quantity: '2.5', minimumQuantity: '1', lotNumber: 'ABC-42', expiryDate: '2026-06-30', notes: 'Top shelf', isFavorite: true }`
   *
   * Returned payload example:
   * - `{ material_id: 12, sector_ids: [7, 8], stock_unit_id: 31, quantity: '2.5', minimum_quantity: '1', lot_number: 'ABC-42', expiry_date: '2026-06-30', notes: 'Top shelf', is_favorite: true }`
   */
  const buildCreateStockPayload = (): CreateInventoryStockPayload | null => {
    const materialId = Number.parseInt(formState.value.materialId, 10)
    const sectorIds = parsePositiveIntegerList(formState.value.sectorIds)
    const stockUnitId = Number.parseInt(formState.value.stockUnitId, 10)
    const quantityValue = parseDecimal(formState.value.quantity)
    const minimumQuantityValue = parseInteger(formState.value.minimumQuantity)

    if (!Number.isInteger(materialId) || materialId <= 0) return null
    if (sectorIds.length === 0) return null
    if (!Number.isInteger(stockUnitId) || stockUnitId <= 0) return null
    if (quantityValue === null || quantityValue <= 0) return null
    if (minimumQuantityValue === null || minimumQuantityValue < 0) return null

    return {
      material_id: materialId,
      sector_ids: sectorIds,
      stock_unit_id: stockUnitId,
      quantity: String(quantityValue),
      minimum_quantity: String(minimumQuantityValue),
      lot_number: formState.value.lotNumber.trim() || null,
      expiry_date: formState.value.expiryDate.trim() || null,
      notes: formState.value.notes.trim() || null,
      is_favorite: formState.value.isFavorite,
    }
  }

  /**
   * Builds the optional material metadata update payload for reagent items.
   *
   * Returned payload examples:
   * - `{ storage_temperature: '4°C' }`
   * - `{ storage_temperature: 'RT', safety_data_sheet: File('sds.pdf') }`
   * - `null`
   */
  const buildReagentMaterialUpdatePayload = (): UpdateInventoryMaterialPayload | null => {
    if (!isSelectedMaterialReagent.value || selectedMaterialId.value <= 0) {
      return null
    }

    const materialDetail = selectedMaterialQuery.data.value
    const nextStorageTemperature = formState.value.reagentStorageTemperature.trim()
    const nextSafetyDataSheet = formState.value.reagentSafetyDataSheet

    const payload: UpdateInventoryMaterialPayload = {}

    if (nextStorageTemperature !== '' && nextStorageTemperature !== (materialDetail?.storage_temperature || '')) {
      payload.storage_temperature = nextStorageTemperature
    }

    if (nextSafetyDataSheet) {
      payload.safety_data_sheet = nextSafetyDataSheet
    }

    return Object.keys(payload).length > 0 ? payload : null
  }

  const submitForm = async (): Promise<void> => {
    if (isSaveDisabled.value) {
      return
    }

    const payload = buildCreateStockPayload()
    if (!payload) {
      toast.add({
        title: 'Please complete required fields',
        color: 'warning',
        duration: 2500,
      })
      return
    }

    try {
      const reagentMaterialPayload = buildReagentMaterialUpdatePayload()

      if (reagentMaterialPayload && selectedMaterialId.value > 0) {
        await inventoryMaterialStore.updateMaterial(selectedMaterialId.value, reagentMaterialPayload)
        await queryClient.invalidateQueries({ queryKey: INVENTORY_MATERIALS_QUERY_KEY })
        await queryClient.invalidateQueries({ queryKey: getInventoryMaterialQueryKey(selectedMaterialId.value) })
      }

      await inventoryStockStore.createStock(payload)
      await queryClient.invalidateQueries({ queryKey: INVENTORY_STOCKS_QUERY_KEY })
      await queryClient.invalidateQueries({ queryKey: INVENTORY_STOCK_PAGES_QUERY_KEY })

      toast.add({
        title: 'Item added',
        color: 'success',
        duration: 2500,
      })

      onSaved()
    } catch (err: unknown) {
      toast.add({
        title: 'Failed to add item',
        description: getErrorMessage(err),
        color: 'error',
        duration: 4000,
      })
    }
  }

  const setReagentSafetyDataSheet = (file: File | null): void => {
    formState.value.reagentSafetyDataSheet = file
  }

  const requestOrderPrefill = (): void => {
    if (isOrderPrefillDisabled.value || isSubmitting.value) {
      return
    }

    void (async (): Promise<void> => {
      const order = selectedOrder.value
      if (!order) {
        toast.add({
          title: 'Select a valid order first',
          color: 'warning',
          duration: 2500,
        })
        return
      }

      formState.value.materialId = String(order.material.id)
      formState.value.stockUnitId = ''

      await nextTick()

      const materialRefetch = await selectedMaterialQuery.refetch()
      const materialDetail = materialRefetch.data ?? selectedMaterialQuery.data.value

      if (!materialDetail || materialDetail.id !== order.material.id) {
        toast.add({
          title: 'Failed to load material units for selected order',
          color: 'error',
          duration: 3500,
        })
        return
      }

      const stockUnits = (materialDetail.units ?? []).filter((unit) => unit.is_stock_unit)
      if (stockUnits.length === 0) {
        toast.add({
          title: 'Selected material has no stock unit',
          color: 'warning',
          duration: 3500,
        })
        return
      }

      const stockUnit = stockUnits[0]
      formState.value.stockUnitId = String(stockUnit.id)

      const orderAmount = parseDecimal(order.amount)
      const orderBaseUnitsPerUnit = parseDecimal(order.order_unit.base_units_per_unit)
      const stockBaseUnitsPerUnit = parseDecimal(stockUnit.base_units_per_unit)

      const canConvertQuantity =
        orderAmount !== null &&
        orderBaseUnitsPerUnit !== null &&
        stockBaseUnitsPerUnit !== null &&
        orderBaseUnitsPerUnit > 0 &&
        stockBaseUnitsPerUnit > 0

      if (!canConvertQuantity) {
        formState.value.quantity = ''
        toast.add({
          title: 'Could not convert order amount to stock unit',
          color: 'warning',
          duration: 3500,
        })
      } else {
        const quantityInBaseUnits = orderAmount * orderBaseUnitsPerUnit
        const quantityInStockUnits = quantityInBaseUnits / stockBaseUnitsPerUnit
        formState.value.quantity = formatDecimal(quantityInStockUnits)
      }

      if (formState.value.notes.trim() === '') {
        const orderNotes = order.notes?.trim() || ''
        formState.value.notes = orderNotes === '' ? `Order #${order.id}` : `Order #${order.id}: ${orderNotes}`
      }

      toast.add({
        title: 'Form prefilled from order',
        color: 'success',
        duration: 2500,
      })
    })()
  }

  return {
    formState,
    selectedOrderId,
    selectedItemTypeId,
    selectedMaterialId,
    selectedRoomId,
    roomOptions,
    supplierDetails,
    filteredMaterials,
    filteredOrders,
    sortedMaterials,
    sortedOrders,
    orderPrefillOptions,
    sectorOptions,
    stockUnitOptions,
    isStockUnitsLoading,
    materialsQuery,
    lookupsQuery,
    ordersQuery,
    materialsErrorMessage,
    sectorsErrorMessage,
    stockUnitsErrorMessage,
    ordersErrorMessage,
    validationMessages,
    isFormValid,
    isOrderPrefillDisabled,
    isSaveDisabled,
    isSubmitting,
    isSelectedMaterialReagent,
    requestOrderPrefill,
    setReagentSafetyDataSheet,
    submitForm,
  }
}
