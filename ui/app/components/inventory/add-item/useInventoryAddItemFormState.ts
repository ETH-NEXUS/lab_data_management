import { computed, ref, watch, type ComputedRef, type Ref } from 'vue'
import { useInventoryLookupsQuery } from '~/composables/inventory/useInventoryLookupQuery'
import {
  useInventoryMaterialQuery,
  useInventoryMaterialsQuery,
} from '~/composables/inventory/useInventoryMaterialQuery'
import { useInventoryOrdersQuery } from '~/composables/inventory/useInventoryOrderQuery'
import type { InventoryMaterialListItem, InventoryOrderListItem } from '~/types/inventory'
import { buildInitialFormState } from '~/components/inventory/add-item/inventoryAddItemForm.utils'
import type { AddItemFormState } from '~/components/inventory/add-item/inventoryAddItemForm.utils'
import { getErrorMessage } from '~/utils/errors'

type UseInventoryAddItemFormStateParams = {
  open: ComputedRef<boolean>
}

type UseInventoryAddItemFormStateResult = {
  formState: Ref<AddItemFormState>
  selectedOrderId: Ref<string>
  selectedItemTypeId: ComputedRef<number>
  selectedMaterialId: ComputedRef<number>
  selectedRoomId: ComputedRef<number>
  roomOptions: ComputedRef<Array<{ id: number; label: string }>>
  selectedMaterialQuery: ReturnType<typeof useInventoryMaterialQuery>
  sortedMaterials: ComputedRef<InventoryMaterialListItem[]>
  filteredMaterials: ComputedRef<InventoryMaterialListItem[]>
  sortedOrders: ComputedRef<InventoryOrderListItem[]>
  filteredOrders: ComputedRef<InventoryOrderListItem[]>
  selectedOrder: ComputedRef<InventoryOrderListItem | null>
  orderPrefillOptions: ComputedRef<Array<{ id: number; label: string }>>
  sectorOptions: ComputedRef<Array<{ id: number; label: string }>>
  brandOptions: ComputedRef<Array<{ id: number; label: string }>>
  manufacturerOptions: ComputedRef<Array<{ id: number; label: string }>>
  vendorOptions: ComputedRef<Array<{ id: number; label: string }>>
  stockUnitOptions: ComputedRef<Array<{ id: number; label: string }>>
  isStockUnitsLoading: ComputedRef<boolean>
  materialsQuery: ReturnType<typeof useInventoryMaterialsQuery>
  lookupsQuery: ReturnType<typeof useInventoryLookupsQuery>
  ordersQuery: ReturnType<typeof useInventoryOrdersQuery>
  materialsErrorMessage: ComputedRef<string | null>
  sectorsErrorMessage: ComputedRef<string | null>
  stockUnitsErrorMessage: ComputedRef<string | null>
  ordersErrorMessage: ComputedRef<string | null>
}

/**
 * Owns the reactive add-item draft state and derived query-backed option lists.
 */
export const useInventoryAddItemFormState = ({
  open,
}: UseInventoryAddItemFormStateParams): UseInventoryAddItemFormStateResult => {
  const materialsQuery = useInventoryMaterialsQuery()
  const lookupsQuery = useInventoryLookupsQuery()
  const ordersQuery = useInventoryOrdersQuery()

  const formState = ref<AddItemFormState>(buildInitialFormState())
  const selectedOrderId = ref('')

  /**
   * Resolves the currently selected item type identifier from the draft form.
   *
   * Input examples:
   * - `{ itemTypeId: '4' }`
   * - `{ itemTypeId: '' }`
   *
   * Returned examples:
   * - `4`
   * - `0`
   */
  const selectedItemTypeId = computed<number>(() => {
    const parsedId = Number.parseInt(formState.value.itemTypeId, 10)
    return Number.isInteger(parsedId) && parsedId > 0 ? parsedId : 0
  })

  const selectedMaterialId = computed<number>(() => {
    const parsedId = Number.parseInt(formState.value.materialId, 10)
    return Number.isInteger(parsedId) && parsedId > 0 ? parsedId : 0
  })

  /**
   * Resolves the selected room identifier from the location step.
   *
   * Input examples:
   * - `{ roomId: '3' }`
   * - `{ roomId: '' }`
   *
   * Returned examples:
   * - `3`
   * - `0`
   */
  const selectedRoomId = computed<number>(() => {
    const parsedId = Number.parseInt(formState.value.roomId, 10)
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

  /**
   * Narrows materials to the selected item type for the first step of the flow.
   *
   * Input examples:
   * - `selectedItemTypeId = 0`
   * - `selectedItemTypeId = 5`
   *
   * Returned examples:
   * - `[]`
   * - `[{ id: 12, item_type: { id: 5, label: 'plate' } }]`
   */
  const filteredMaterials = computed<InventoryMaterialListItem[]>(() => {
    if (selectedItemTypeId.value <= 0) {
      return []
    }

    return sortedMaterials.value.filter((material) => {
      return material.item_type?.id === selectedItemTypeId.value
    })
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

  /**
   * Narrows order-prefill candidates to orders whose material belongs to the selected item type.
   *
   * Input examples:
   * - `selectedItemTypeId = 0`
   * - `selectedItemTypeId = 9`
   *
   * Returned examples:
   * - `[]`
   * - `[{ id: 41, material: { item_type: { id: 9, label: 'bottle' } } }]`
   */
  const filteredOrders = computed<InventoryOrderListItem[]>(() => {
    if (selectedItemTypeId.value <= 0) {
      return []
    }

    return sortedOrders.value.filter((order) => {
      return order.material.item_type?.id === selectedItemTypeId.value
    })
  })

  const selectedOrder = computed<InventoryOrderListItem | null>(() => {
    const orderId = Number.parseInt(selectedOrderId.value, 10)
    if (!Number.isInteger(orderId) || orderId <= 0) {
      return null
    }

    return filteredOrders.value.find((order) => order.id === orderId) ?? null
  })

  const orderPrefillOptions = computed<Array<{ id: number; label: string }>>(() => {
    return filteredOrders.value.map((order) => {
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

  const roomOptions = computed<Array<{ id: number; label: string }>>(() => {
    const rooms = [...(lookupsQuery.data.value?.rooms ?? [])]

    rooms.sort((leftRoom, rightRoom) => {
      const leftLabel = (leftRoom.label || leftRoom.name || '').toLowerCase()
      const rightLabel = (rightRoom.label || rightRoom.name || '').toLowerCase()
      return leftLabel.localeCompare(rightLabel)
    })

    return rooms.map((room) => ({
      id: room.id,
      label: room.label || room.name || `Room #${room.id}`,
    }))
  })

  // The room is already chosen via a separate control before this list is shown,
  // so the option label only needs the sector's own name, not "Room / Sector".
  const sectorOptions = computed<Array<{ id: number; label: string }>>(() => {
    const sectorsInSelectedRoom = (lookupsQuery.data.value?.sectors ?? []).filter(
      (sector) => selectedRoomId.value > 0 && sector.room?.id === selectedRoomId.value,
    )

    sectorsInSelectedRoom.sort((leftSector, rightSector) => {
      const leftLabel = (leftSector.name || '').toLowerCase()
      const rightLabel = (rightSector.name || '').toLowerCase()
      return leftLabel.localeCompare(rightLabel)
    })

    return sectorsInSelectedRoom.map((sector) => ({
      id: sector.id,
      label: sector.name || `Sector #${sector.id}`,
    }))
  })

  // Batch 1: brand options for the general "Additional material details" section.
  const brandOptions = computed<Array<{ id: number; label: string }>>(() => {
    const brands = [...(lookupsQuery.data.value?.brands ?? [])]

    brands.sort((leftBrand, rightBrand) => {
      const leftLabel = (leftBrand.label || leftBrand.name || '').toLowerCase()
      const rightLabel = (rightBrand.label || rightBrand.name || '').toLowerCase()
      return leftLabel.localeCompare(rightLabel)
    })

    return brands.map((brand) => ({
      id: brand.id,
      label: brand.label || brand.name || `Brand #${brand.id}`,
    }))
  })

  // Batch 2: manufacturer/vendor options for the "Additional material details" section.
  const manufacturerOptions = computed<Array<{ id: number; label: string }>>(() => {
    const manufacturers = [...(lookupsQuery.data.value?.manufacturers ?? [])]

    manufacturers.sort((leftManufacturer, rightManufacturer) => {
      const leftLabel = (leftManufacturer.label || leftManufacturer.name || '').toLowerCase()
      const rightLabel = (rightManufacturer.label || rightManufacturer.name || '').toLowerCase()
      return leftLabel.localeCompare(rightLabel)
    })

    return manufacturers.map((manufacturer) => ({
      id: manufacturer.id,
      label: manufacturer.label || manufacturer.name || `Manufacturer #${manufacturer.id}`,
    }))
  })

  const vendorOptions = computed<Array<{ id: number; label: string }>>(() => {
    const vendors = [...(lookupsQuery.data.value?.vendors ?? [])]

    vendors.sort((leftVendor, rightVendor) => {
      const leftLabel = (leftVendor.label || leftVendor.name || '').toLowerCase()
      const rightLabel = (rightVendor.label || rightVendor.name || '').toLowerCase()
      return leftLabel.localeCompare(rightLabel)
    })

    return vendors.map((vendor) => ({
      id: vendor.id,
      label: vendor.label || vendor.name || `Vendor #${vendor.id}`,
    }))
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

  watch(open, (isOpen) => {
    if (!isOpen) {
      return
    }

    formState.value = buildInitialFormState()
    selectedOrderId.value = ''
  })

  watch(
    () => formState.value.itemTypeId,
    () => {
      formState.value.materialId = ''
      formState.value.stockUnitId = ''
      selectedOrderId.value = ''
    },
  )

  watch(
    () => formState.value.materialId,
    () => {
      formState.value.reagentStorageTemperature = ''
      formState.value.reagentSafetyDataSheet = null
      formState.value.additionalBrandId = ''
      formState.value.additionalDefaultCost = ''
      formState.value.additionalManufacturerId = ''
      formState.value.additionalVendorId = ''
      formState.value.additionalManufacturerCatalogNumber = ''
      formState.value.additionalVendorCatalogNumber = ''
      formState.value.additionalCapacityValue = ''
      formState.value.additionalCapacityUnit = ''
      formState.value.additionalDescription = ''
      formState.value.additionalSerialNumber = ''
      formState.value.additionalOrderNumber = ''
      formState.value.additionalLifetimeDays = ''
      formState.value.additionalIsActive = ''
      formState.value.stockUnitId = ''
    },
  )

  watch(
    () => selectedMaterialQuery.data.value,
    (materialDetail) => {
      if (!materialDetail || materialDetail.id !== selectedMaterialId.value) {
        return
      }

      formState.value.reagentStorageTemperature = materialDetail.storage_temperature || ''
      formState.value.reagentSafetyDataSheet = null
      formState.value.additionalBrandId = materialDetail.brand?.id ? String(materialDetail.brand.id) : ''
      formState.value.additionalDefaultCost = materialDetail.default_cost || ''
      formState.value.additionalManufacturerId = materialDetail.manufacturer?.id
        ? String(materialDetail.manufacturer.id)
        : ''
      formState.value.additionalVendorId = materialDetail.vendor?.id ? String(materialDetail.vendor.id) : ''
      formState.value.additionalManufacturerCatalogNumber = materialDetail.manufacturer_catalog_number || ''
      formState.value.additionalVendorCatalogNumber = materialDetail.vendor_catalog_number || ''
      formState.value.additionalCapacityValue = materialDetail.capacity_value || ''
      formState.value.additionalCapacityUnit = materialDetail.capacity_unit || ''
      formState.value.additionalDescription = materialDetail.description || ''
      formState.value.additionalSerialNumber = materialDetail.serial_number || ''
      formState.value.additionalOrderNumber = materialDetail.order_number || ''
      formState.value.additionalLifetimeDays =
        materialDetail.lifetime_days == null ? '' : String(materialDetail.lifetime_days)
      // additionalIsActive intentionally stays '' (not prefilled): '' means
      // "don't change" for this tri-state field, unlike the other text fields above.
    },
  )

  watch(
    () => formState.value.roomId,
    () => {
      formState.value.sectorIds = []
    },
  )

  return {
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
  }
}
