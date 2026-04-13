import { computed, nextTick, ref, watch, type ComputedRef } from 'vue'
import { useQueryClient } from '@tanstack/vue-query'
import { useInventoryLookupsQuery } from '~/composables/inventory/useInventoryLookupQuery'
import {
  useInventoryMaterialQuery,
  useInventoryMaterialsQuery,
} from '~/composables/inventory/useInventoryMaterialQuery'
import { useInventoryOrdersQuery } from '~/composables/inventory/useInventoryOrderQuery'
import { useInventoryStockStore } from '~/stores/inventory/InventoryStock'
import {
  INVENTORY_STOCKS_QUERY_KEY,
  type CreateInventoryStockPayload,
  type InventoryMaterialListItem,
  type InventoryOrderListItem,
} from '~/types/inventory'
import { getErrorMessage } from '~/utils/errors'

export type AddItemFormState = {
  materialId: string
  sectorId: string
  stockUnitId: string
  quantity: string
  minimumQuantity: string
  lotNumber: string
  expiryDate: string
  notes: string
  isFavorite: boolean
}

type UseInventoryAddItemFormParams = {
  open: ComputedRef<boolean>
  onSaved: () => void
}

/**
 * Builds empty draft values for the add-item stock form.
 *
 * Returned data example:
 * - `{ materialId: '', sectorId: '', stockUnitId: '', quantity: '', minimumQuantity: '', lotNumber: '', expiryDate: '', notes: '', isFavorite: false }`
 */
const buildInitialFormState = (): AddItemFormState => ({
  materialId: '',
  sectorId: '',
  stockUnitId: '',
  quantity: '',
  minimumQuantity: '',
  lotNumber: '',
  expiryDate: '',
  notes: '',
  isFavorite: false,
})

/**
 * Parses decimal-like text into a finite number.
 *
 * Input examples:
 * - `'2.5'`
 * - `'96'`
 * - `null`
 *
 * Returned examples:
 * - `2.5`
 * - `96`
 * - `null`
 */
const parseDecimal = (value: string | null | undefined): number | null => {
  if (value == null) {
    return null
  }

  const parsedValue = Number.parseFloat(value.trim())
  if (!Number.isFinite(parsedValue)) {
    return null
  }

  return parsedValue
}

/**
 * Formats decimal values for API-compatible string fields.
 *
 * Input examples:
 * - `2`
 * - `2.041666666`
 *
 * Returned examples:
 * - `'2'`
 * - `'2.041667'`
 */
const formatDecimal = (value: number): string => {
  const fixedValue = value.toFixed(6)
  const normalizedValue = fixedValue.replace(/\.?0+$/, '')
  return normalizedValue === '' ? '0' : normalizedValue
}

/**
 * Manages add-item modal form state, dependent lookup data, validation, and stock creation submit.
 */
export const useInventoryAddItemForm = ({ open, onSaved }: UseInventoryAddItemFormParams) => {
  const toast = useToast()
  const queryClient = useQueryClient()
  const inventoryStockStore = useInventoryStockStore()
  const materialsQuery = useInventoryMaterialsQuery()
  const lookupsQuery = useInventoryLookupsQuery()
  const ordersQuery = useInventoryOrdersQuery()
  const isSubmitting = computed<boolean>(() => inventoryStockStore.isCreatingStock)

  const formState = ref<AddItemFormState>(buildInitialFormState())
  const selectedOrderId = ref('')

  const selectedMaterialId = computed<number>(() => {
    const parsedId = Number.parseInt(formState.value.materialId, 10)
    return Number.isInteger(parsedId) && parsedId > 0 ? parsedId : 0
  })

  const selectedMaterialQuery = useInventoryMaterialQuery(selectedMaterialId)

  const sortedMaterials = computed<InventoryMaterialListItem[]>(() => {
    const materials = [...(materialsQuery.data.value ?? [])]
    materials.sort((leftMaterial, rightMaterial) => {
      const leftLabel = (leftMaterial.label || leftMaterial.product_name || '').toLowerCase()
      const rightLabel = (rightMaterial.label || rightMaterial.product_name || '').toLowerCase()
      return leftLabel.localeCompare(rightLabel)
    })
    return materials
  })

  const sortedOrders = computed<InventoryOrderListItem[]>(() => {
    const orders = [...(ordersQuery.data.value ?? [])]
    orders.sort((leftOrder, rightOrder) => {
      const leftDate = new Date(leftOrder.order_date).getTime()
      const rightDate = new Date(rightOrder.order_date).getTime()
      return rightDate - leftDate
    })
    return orders
  })

  const selectedOrder = computed<InventoryOrderListItem | null>(() => {
    const orderId = Number.parseInt(selectedOrderId.value, 10)
    if (!Number.isInteger(orderId) || orderId <= 0) {
      return null
    }

    return sortedOrders.value.find((order) => order.id === orderId) ?? null
  })

  const orderPrefillOptions = computed<Array<{ id: number; label: string }>>(() => {
    return sortedOrders.value.map((order) => {
      const materialLabel = order.material.label || order.material.product_name || `Material #${order.material.id}`
      const unitLabel =
        order.order_unit.display_name || order.order_unit.unit?.label || order.order_unit.unit?.name || ''
      const statusLabel = order.status_label || order.status
      return {
        id: order.id,
        label: `#${order.id} · ${materialLabel} · ${order.amount} ${unitLabel}`.trim() + ` · ${statusLabel}`,
      }
    })
  })

  const sectorOptions = computed<Array<{ id: number; label: string }>>(() => {
    const sectors = [...(lookupsQuery.data.value?.sectors ?? [])]
    sectors.sort((leftSector, rightSector) => {
      const leftRoomLabel = (leftSector.room?.label || leftSector.room?.name || '').toLowerCase()
      const rightRoomLabel = (rightSector.room?.label || rightSector.room?.name || '').toLowerCase()
      if (leftRoomLabel !== rightRoomLabel) {
        return leftRoomLabel.localeCompare(rightRoomLabel)
      }

      const leftSectorLabel = (leftSector.label || leftSector.name || '').toLowerCase()
      const rightSectorLabel = (rightSector.label || rightSector.name || '').toLowerCase()
      return leftSectorLabel.localeCompare(rightSectorLabel)
    })

    return sectors.map((sector) => {
      const roomLabel = sector.room?.label || sector.room?.name || 'Unknown room'
      const sectorLabel = sector.name || `Sector #${sector.id}`
      return {
        id: sector.id,
        label: sector.label || `${roomLabel} / ${sectorLabel}`,
      }
    })
  })

  const stockUnitOptions = computed<Array<{ id: number; label: string }>>(() => {
    const material = selectedMaterialQuery.data.value
    const units = material?.units ?? []

    return units
      .filter((unit) => unit.is_stock_unit)
      .map((unit) => ({
        id: unit.id,
        label: unit.display_name || unit.unit?.label || unit.unit?.name || `Unit #${unit.id}`,
      }))
  })

  const isStockUnitsLoading = computed<boolean>(() => {
    return selectedMaterialId.value > 0 && selectedMaterialQuery.isPending.value
  })

  const materialsErrorMessage = computed<string | null>(() => {
    const err = materialsQuery.error.value
    if (!err) {
      return null
    }
    return getErrorMessage(err)
  })

  const sectorsErrorMessage = computed<string | null>(() => {
    const err = lookupsQuery.error.value
    if (!err) {
      return null
    }
    return getErrorMessage(err)
  })

  const stockUnitsErrorMessage = computed<string | null>(() => {
    const err = selectedMaterialQuery.error.value
    if (!err) {
      return null
    }
    return getErrorMessage(err)
  })

  const ordersErrorMessage = computed<string | null>(() => {
    const err = ordersQuery.error.value
    if (!err) {
      return null
    }
    return getErrorMessage(err)
  })

  const validationMessages = computed<string[]>(() => {
    const messages: string[] = []
    const materialId = Number.parseInt(formState.value.materialId, 10)
    const sectorId = Number.parseInt(formState.value.sectorId, 10)
    const stockUnitId = Number.parseInt(formState.value.stockUnitId, 10)
    const quantityValue = Number.parseFloat(formState.value.quantity.trim())
    const minimumQuantityValue = Number.parseFloat(formState.value.minimumQuantity.trim())

    if (!Number.isInteger(materialId) || materialId <= 0) {
      messages.push('Material is required.')
    }
    if (!Number.isInteger(sectorId) || sectorId <= 0) {
      messages.push('Sector is required.')
    }
    if (!Number.isInteger(stockUnitId) || stockUnitId <= 0) {
      messages.push('Stock unit is required.')
    }
    if (!Number.isFinite(quantityValue) || quantityValue <= 0) {
      messages.push('Quantity must be greater than zero.')
    }
    if (!Number.isFinite(minimumQuantityValue) || minimumQuantityValue < 0) {
      messages.push('Minimum quantity must be zero or greater.')
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

  watch(open, (isOpen) => {
    if (!isOpen) {
      return
    }

    formState.value = buildInitialFormState()
    selectedOrderId.value = ''
  })

  watch(
    () => formState.value.materialId,
    () => {
      formState.value.stockUnitId = ''
    },
  )

  /**
   * Builds API payload for inventory stock creation from form draft state.
   *
   * Accepted draft example:
   * - `{ materialId: '12', sectorId: '7', stockUnitId: '31', quantity: '2.5', minimumQuantity: '1', lotNumber: 'ABC-42', expiryDate: '2026-06-30', notes: 'Top shelf', isFavorite: true }`
   *
   * Returned payload example:
   * - `{ material_id: 12, sector_id: 7, stock_unit_id: 31, quantity: '2.5', minimum_quantity: '1', lot_number: 'ABC-42', expiry_date: '2026-06-30', notes: 'Top shelf', is_favorite: true }`
   */
  const buildCreateStockPayload = (): CreateInventoryStockPayload | null => {
    const materialId = Number.parseInt(formState.value.materialId, 10)
    const sectorId = Number.parseInt(formState.value.sectorId, 10)
    const stockUnitId = Number.parseInt(formState.value.stockUnitId, 10)
    const quantityValue = Number.parseFloat(formState.value.quantity.trim())
    const minimumQuantityValue = Number.parseFloat(formState.value.minimumQuantity.trim())

    if (!Number.isInteger(materialId) || materialId <= 0) return null
    if (!Number.isInteger(sectorId) || sectorId <= 0) return null
    if (!Number.isInteger(stockUnitId) || stockUnitId <= 0) return null
    if (!Number.isFinite(quantityValue) || quantityValue <= 0) return null
    if (!Number.isFinite(minimumQuantityValue) || minimumQuantityValue < 0) return null

    return {
      material_id: materialId,
      sector_id: sectorId,
      stock_unit_id: stockUnitId,
      quantity: String(quantityValue),
      minimum_quantity: String(minimumQuantityValue),
      lot_number: formState.value.lotNumber.trim() || null,
      expiry_date: formState.value.expiryDate.trim() || null,
      notes: formState.value.notes.trim() || null,
      is_favorite: formState.value.isFavorite,
    }
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
      await inventoryStockStore.createStock(payload)
      await queryClient.invalidateQueries({ queryKey: INVENTORY_STOCKS_QUERY_KEY })

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
    selectedMaterialId,
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
    requestOrderPrefill,
    submitForm,
  }
}
