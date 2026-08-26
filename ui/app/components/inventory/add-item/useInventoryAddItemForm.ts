import { computed, nextTick, ref, watch, type ComputedRef } from 'vue'
import { useQueryClient } from '@tanstack/vue-query'
import { useInventoryMaterialStore } from '~/stores/inventory/InventoryMaterialStore'
import { useInventoryStockStore } from '~/stores/inventory/InventoryStock'
import { useInventoryStockTablePreferenceStore } from '~/stores/inventory/InventoryStockTablePreferenceStore'
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
  const stockTablePreferenceStore = useInventoryStockTablePreferenceStore()
  const isSubmitting = computed<boolean>(() => {
    return inventoryStockStore.isCreatingStock || inventoryMaterialStore.isUpdatingMaterial
  })
  const appliedSourceOrderId = ref<number | null>(null)
  const appliedSourceMaterialId = ref<number | null>(null)

  const resetAppliedSourceOrder = (): void => {
    appliedSourceOrderId.value = null
    appliedSourceMaterialId.value = null
  }

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
    brandOptions,
    manufacturerOptions,
    vendorOptions,
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
    // Optional: left blank, it defaults to 0. Only flag it if something was typed and it's invalid.
    if (formState.value.minimumQuantity.trim() !== '' && (minimumQuantityValue === null || minimumQuantityValue < 0)) {
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
   * - `{ materialId: '12', sectorIds: ['7', '8'], stockUnitId: '31', quantity: '2.5', minimumQuantity: '1', lotNumber: 'ABC-42', expiryDate: '2026-06-30', notes: 'Top shelf' }`
   *
   * Returned payload example:
   * - `{ material_id: 12, sector_ids: [7, 8], stock_unit_id: 31, quantity: '2.5', minimum_quantity: '1', lot_number: 'ABC-42', expiry_date: '2026-06-30', notes: 'Top shelf' }`
   */
  const buildCreateStockPayload = (): CreateInventoryStockPayload | null => {
    const materialId = Number.parseInt(formState.value.materialId, 10)
    const sectorIds = parsePositiveIntegerList(formState.value.sectorIds)
    const stockUnitId = Number.parseInt(formState.value.stockUnitId, 10)
    const quantityValue = parseDecimal(formState.value.quantity)
    // Optional: blank means "use the default of 0", same as the model's own default.
    const isMinimumQuantityBlank = formState.value.minimumQuantity.trim() === ''
    const minimumQuantityValue = isMinimumQuantityBlank ? 0 : parseInteger(formState.value.minimumQuantity)

    if (!Number.isInteger(materialId) || materialId <= 0) return null
    if (sectorIds.length === 0) return null
    if (!Number.isInteger(stockUnitId) || stockUnitId <= 0) return null
    if (quantityValue === null || quantityValue <= 0) return null
    if (minimumQuantityValue === null || minimumQuantityValue < 0) return null

    const shouldLinkSourceOrder =
      appliedSourceOrderId.value !== null &&
      appliedSourceMaterialId.value !== null &&
      appliedSourceMaterialId.value === materialId

    return {
      material_id: materialId,
      sector_ids: sectorIds,
      stock_unit_id: stockUnitId,
      quantity: String(quantityValue),
      minimum_quantity: String(minimumQuantityValue),
      source_order_id: shouldLinkSourceOrder ? appliedSourceOrderId.value : undefined,
      lot_number: formState.value.lotNumber.trim() || null,
      expiry_date: formState.value.expiryDate.trim() || null,
      notes: formState.value.notes.trim() || null,
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

  /**
   * Builds the optional material metadata update payload for the general
   * "Additional material details" section. Every field here is optional:
   * left blank, it's simply skipped and the material is left unchanged.
   *
   * Returned payload examples:
   * - `{ brand_id: 4 }`
   * - `{ default_cost: '12.50' }`
   * - `null`
   */
  const buildAdditionalMaterialUpdatePayload = (): UpdateInventoryMaterialPayload | null => {
    if (selectedMaterialId.value <= 0) {
      return null
    }

    const materialDetail = selectedMaterialQuery.data.value
    const payload: UpdateInventoryMaterialPayload = {}

    const nextBrandId = Number.parseInt(formState.value.additionalBrandId, 10)
    const currentBrandId = materialDetail?.brand?.id ?? null
    if (Number.isInteger(nextBrandId) && nextBrandId > 0 && nextBrandId !== currentBrandId) {
      payload.brand_id = nextBrandId
    }

    const nextDefaultCost = formState.value.additionalDefaultCost.trim()
    if (nextDefaultCost !== '' && nextDefaultCost !== (materialDetail?.default_cost || '')) {
      payload.default_cost = nextDefaultCost
    }

    const nextManufacturerId = Number.parseInt(formState.value.additionalManufacturerId, 10)
    const currentManufacturerId = materialDetail?.manufacturer?.id ?? null
    if (Number.isInteger(nextManufacturerId) && nextManufacturerId > 0 && nextManufacturerId !== currentManufacturerId) {
      payload.manufacturer_id = nextManufacturerId
    }

    const nextVendorId = Number.parseInt(formState.value.additionalVendorId, 10)
    const currentVendorId = materialDetail?.vendor?.id ?? null
    if (Number.isInteger(nextVendorId) && nextVendorId > 0 && nextVendorId !== currentVendorId) {
      payload.vendor_id = nextVendorId
    }

    const nextManufacturerCatalogNumber = formState.value.additionalManufacturerCatalogNumber.trim()
    if (
      nextManufacturerCatalogNumber !== '' &&
      nextManufacturerCatalogNumber !== (materialDetail?.manufacturer_catalog_number || '')
    ) {
      payload.manufacturer_catalog_number = nextManufacturerCatalogNumber
    }

    const nextVendorCatalogNumber = formState.value.additionalVendorCatalogNumber.trim()
    if (nextVendorCatalogNumber !== '' && nextVendorCatalogNumber !== (materialDetail?.vendor_catalog_number || '')) {
      payload.vendor_catalog_number = nextVendorCatalogNumber
    }

    const nextCapacityValue = formState.value.additionalCapacityValue.trim()
    if (nextCapacityValue !== '' && nextCapacityValue !== (materialDetail?.capacity_value || '')) {
      payload.capacity_value = nextCapacityValue
    }

    const nextCapacityUnit = formState.value.additionalCapacityUnit.trim()
    if (nextCapacityUnit !== '' && nextCapacityUnit !== (materialDetail?.capacity_unit || '')) {
      payload.capacity_unit = nextCapacityUnit
    }

    const nextDescription = formState.value.additionalDescription.trim()
    if (nextDescription !== '' && nextDescription !== (materialDetail?.description || '')) {
      payload.description = nextDescription
    }

    const nextSerialNumber = formState.value.additionalSerialNumber.trim()
    if (nextSerialNumber !== '' && nextSerialNumber !== (materialDetail?.serial_number || '')) {
      payload.serial_number = nextSerialNumber
    }

    const nextOrderNumber = formState.value.additionalOrderNumber.trim()
    if (nextOrderNumber !== '' && nextOrderNumber !== (materialDetail?.order_number || '')) {
      payload.order_number = nextOrderNumber
    }

    const nextLifetimeDaysText = formState.value.additionalLifetimeDays.trim()
    if (nextLifetimeDaysText !== '') {
      const nextLifetimeDays = Number.parseInt(nextLifetimeDaysText, 10)
      if (
        Number.isInteger(nextLifetimeDays) &&
        nextLifetimeDays >= 0 &&
        nextLifetimeDays !== materialDetail?.lifetime_days
      ) {
        payload.lifetime_days = nextLifetimeDays
      }
    }

    const nextIsActiveText = formState.value.additionalIsActive
    if (nextIsActiveText !== '') {
      const nextIsActive = nextIsActiveText === 'true'
      if (nextIsActive !== materialDetail?.is_active) {
        payload.is_active = nextIsActive
      }
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
      // Two independent optional payloads (reagent fields, general additional fields)
      // are merged into one PATCH so the material is only ever updated once.
      const reagentMaterialPayload = buildReagentMaterialUpdatePayload()
      const additionalMaterialPayload = buildAdditionalMaterialUpdatePayload()
      const materialUpdatePayload =
        reagentMaterialPayload || additionalMaterialPayload
          ? { ...additionalMaterialPayload, ...reagentMaterialPayload }
          : null

      if (materialUpdatePayload && selectedMaterialId.value > 0) {
        await inventoryMaterialStore.updateMaterial(selectedMaterialId.value, materialUpdatePayload)
        await queryClient.invalidateQueries({ queryKey: INVENTORY_MATERIALS_QUERY_KEY })
        await queryClient.invalidateQueries({ queryKey: getInventoryMaterialQueryKey(selectedMaterialId.value) })
      }

      const createdStock = await inventoryStockStore.createStock(payload)
      await queryClient.invalidateQueries({ queryKey: INVENTORY_STOCKS_QUERY_KEY })
      await queryClient.invalidateQueries({ queryKey: INVENTORY_STOCK_PAGES_QUERY_KEY })

      toast.add({
        title: 'Item added',
        color: 'success',
        duration: 2500,
      })

      await stockTablePreferenceStore.updateGlobalFilterState('')
      onSaved()
      await navigateTo({
        path: '/inventory/all',
        query: {
          preset: 'all',
          stock: String(createdStock.id),
        },
      })
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

      appliedSourceOrderId.value = order.id
      appliedSourceMaterialId.value = order.material.id

      toast.add({
        title: 'Form prefilled from order',
        color: 'success',
        duration: 2500,
      })
    })()
  }

  watch(open, (isOpen) => {
    if (!isOpen) {
      return
    }

    resetAppliedSourceOrder()
  })

  watch(
    () => selectedOrderId.value,
    (nextSelectedOrderId) => {
      if (nextSelectedOrderId.trim() !== '') {
        resetAppliedSourceOrder()
      }
    },
  )

  watch(
    () => formState.value.materialId,
    (nextMaterialId) => {
      if (
        appliedSourceMaterialId.value !== null &&
        nextMaterialId.trim() !== '' &&
        nextMaterialId !== String(appliedSourceMaterialId.value)
      ) {
        resetAppliedSourceOrder()
      }
    },
  )

  return {
    formState,
    selectedOrderId,
    selectedItemTypeId,
    selectedMaterialId,
    selectedRoomId,
    roomOptions,
    filteredMaterials,
    filteredOrders,
    sortedMaterials,
    sortedOrders,
    orderPrefillOptions,
    sectorOptions,
    brandOptions,
    manufacturerOptions,
    vendorOptions,
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
