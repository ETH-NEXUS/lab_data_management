import { computed, ref, watch } from 'vue'
import { useI18n } from '#imports'
import type {
  InventoryMaterialDetail,
  InventoryOrderListItem,
  InventoryStockListItem,
  InventoryUsageListItem,
} from '~/types/inventory'
import { formatDateTime } from '~/utils/dateTime'

type DrawerField = {
  label: string
  value: string
  wide?: boolean
}

export type InventoryStockDetailsPanelProps = {
  open: boolean
  stock: InventoryStockListItem | null
  materialDetail: InventoryMaterialDetail | null
  usageEntries: InventoryUsageListItem[]
  orderEntries: InventoryOrderListItem[]
}

/**
 * Builds all computed display fields used by the inventory stock details panel.
 *
 * Accepted props example:
 * - `{ open: true, stock: { id: 14, ... }, materialDetail: { id: 3, ... }, usageEntries: [], orderEntries: [] }`
 */
export const useInventoryStockDetailsPanel = (props: Readonly<InventoryStockDetailsPanelProps>) => {
  const { t } = useI18n()

  const isMetadataExpanded = ref(false)

  watch(
    () => props.open,
    (isOpen) => {
      if (isOpen) {
        isMetadataExpanded.value = false
      }
    },
  )

  /**
   * Converts optional values into a readable fallback-safe string.
   *
   * Input examples:
   * - `'Vendor A'`
   * - `null`
   * - `''`
   *
   * Returned examples:
   * - `'Vendor A'`
   * - `'—'`
   */
  const toDisplayValue = (value: unknown): string => {
    if (value == null) {
      return t('inventory.stock_drawer.values.none')
    }

    const textValue = String(value).trim()
    if (textValue === '') {
      return t('inventory.stock_drawer.values.none')
    }

    return textValue
  }

  /**
   * Normalizes decimal strings for compact UI display.
   *
   * Input examples:
   * - `'10.000000'`
   * - `'3.250000'`
   * - `'5'`
   *
   * Returned examples:
   * - `'10'`
   * - `'3.25'`
   * - `'5'`
   */
  const formatDecimalForDisplay = (value: string | null | undefined): string => {
    const rawValue = value?.trim() || ''
    if (rawValue === '') {
      return ''
    }

    if (!/^-?\d+(\.\d+)?$/.test(rawValue)) {
      return rawValue
    }

    if (!rawValue.includes('.')) {
      return rawValue
    }

    const normalizedValue = rawValue.replace(/\.?0+$/, '')
    return normalizedValue === '-0' ? '0' : normalizedValue
  }

  /**
   * Resolves one location label from stock values.
   *
   * Accepted stock example:
   * - `{ location_label: 'Room A / Shelf 3' }`
   * - `{ location_label: 'C75 / 3.1, 3.2' }`
   *
   * Returned examples:
   * - `'Room A / Shelf 3'`
   * - `'C75 / 3.1, 3.2'`
   * - `'—'`
   */
  const locationLabel = computed<string>(() => {
    const stock = props.stock
    if (!stock) {
      return t('inventory.stock_drawer.values.none')
    }

    if (stock.location_label?.trim()) {
      return stock.location_label.trim()
    }

    const roomLabel =
      stock.room?.label ?? stock.room?.name ?? stock.sector?.room?.label ?? stock.sector?.room?.name ?? ''
    const sectorLabel = stock.sector?.label ?? stock.sector?.name ?? ''

    if (roomLabel !== '' && sectorLabel !== '') {
      return `${roomLabel} / ${sectorLabel}`
    }

    if (roomLabel !== '') {
      return roomLabel
    }

    if (sectorLabel !== '') {
      return sectorLabel
    }

    return t('inventory.stock_drawer.values.none')
  })

  const inventoryStatusLabel = computed<string>(() => {
    if (!props.stock) {
      return t('inventory.stock_drawer.values.none')
    }

    if (props.stock.inventory_status === 'out_of_stock') {
      return t('inventory.stock_table.status_labels.out_of_stock')
    }

    return props.stock.inventory_status === 'low'
      ? t('inventory.stock_table.status_labels.low')
      : t('inventory.stock_table.status_labels.in_stock')
  })

  const inventoryStatusColor = computed<'error' | 'warning' | 'success'>(() => {
    if (props.stock?.inventory_status === 'out_of_stock') {
      return 'error'
    }

    return props.stock?.inventory_status === 'low' ? 'warning' : 'success'
  })

  // Reagent-only field on the material. Surfaced as an emphasized badge (not
  // buried in the material identity grid) because handling it wrong matters.
  const storageTemperatureLabel = computed<string | null>(() => {
    const stock = props.stock
    if (!stock) {
      return null
    }

    const label = stock.material.storage_temperature_label || stock.material.storage_temperature || ''
    return label.trim() === '' ? null : label.trim()
  })

  const stockUnitLabel = computed<string>(() => {
    if (!props.stock) {
      return t('inventory.stock_drawer.values.none')
    }

    return toDisplayValue(
      props.stock.stock_unit.display_name || props.stock.stock_unit.unit?.label || props.stock.stock_unit.unit?.name,
    )
  })

  const quantityLabel = computed<string>(() => {
    if (!props.stock) {
      return t('inventory.stock_drawer.values.none')
    }

    return `${toDisplayValue(formatDecimalForDisplay(props.stock.quantity))} ${stockUnitLabel.value}`.trim()
  })

  const minimumQuantityLabel = computed<string>(() => {
    if (!props.stock) {
      return t('inventory.stock_drawer.values.none')
    }

    return `${toDisplayValue(formatDecimalForDisplay(props.stock.minimum_quantity))} ${stockUnitLabel.value}`.trim()
  })

  const baseUnitLabel = computed<string>(() => {
    const stock = props.stock
    if (!stock) {
      return t('inventory.stock_drawer.values.none')
    }

    // Use material detail units when available.
    const materialUnits = props.materialDetail?.units ?? []
    const baseUnit = materialUnits.find((unit) => unit.is_base_unit)
    if (baseUnit) {
      return toDisplayValue(baseUnit.display_name || baseUnit.unit?.label || baseUnit.unit?.name)
    }

    // Fallback: stock unit can already be base.
    if (stock.stock_unit.is_base_unit) {
      return toDisplayValue(
        stock.stock_unit.display_name || stock.stock_unit.unit?.label || stock.stock_unit.unit?.name,
      )
    }

    // Fallback to stock unit label.
    return toDisplayValue(stock.stock_unit.unit?.label || stock.stock_unit.unit?.name)
  })

  const attributeLabel = computed<string>(() => {
    const stock = props.stock
    if (!stock) {
      return t('inventory.stock_drawer.values.none')
    }

    const labels: string[] = []
    for (const attribute of stock.material.attributes ?? []) {
      const label = attribute.label || attribute.name
      if (label && label.trim() !== '') {
        labels.push(label)
      }
    }

    if (labels.length === 0) {
      return t('inventory.stock_drawer.values.no_attributes')
    }

    return labels.join(', ')
  })

  const capacityWithUnitLabel = computed<string>(() => {
    const stock = props.stock
    if (!stock) {
      return t('inventory.stock_drawer.values.none')
    }

    const capacityValue = toDisplayValue(formatDecimalForDisplay(stock.material.capacity_value))
    const capacityUnit = toDisplayValue(stock.material.capacity_unit)

    if (
      capacityValue === t('inventory.stock_drawer.values.none') &&
      capacityUnit === t('inventory.stock_drawer.values.none')
    ) {
      return t('inventory.stock_drawer.values.none')
    }

    if (capacityValue === t('inventory.stock_drawer.values.none')) {
      return capacityUnit
    }

    if (capacityUnit === t('inventory.stock_drawer.values.none')) {
      return capacityValue
    }

    return `${capacityValue} ${capacityUnit}`.trim()
  })

  const operationalFields = computed<DrawerField[]>(() => {
    const stock = props.stock
    if (!stock) return []

    return [
      {
        label: t('inventory.stock_drawer.fields.product_name'),
        value: toDisplayValue(stock.material.product_name),
        wide: true,
      },
      {
        label: t('inventory.stock_drawer.fields.quantity'),
        value: quantityLabel.value,
      },
      {
        label: t('inventory.stock_drawer.fields.minimum_quantity'),
        value: minimumQuantityLabel.value,
      },
      {
        label: t('inventory.stock_drawer.fields.location'),
        value: locationLabel.value,
        wide: true,
      },
      {
        label: t('inventory.stock_drawer.fields.lot_number'),
        value: toDisplayValue(stock.lot_number),
      },
      {
        label: t('inventory.stock_drawer.fields.expiry_date'),
        value: stock.expiry_date ? stock.expiry_date : t('inventory.stock_drawer.values.none'),
      },
      {
        label: t('inventory.stock_drawer.fields.notes'),
        value: stock.notes || t('inventory.stock_drawer.values.no_notes'),
        wide: true,
      },
    ]
  })

  const unitConversionFields = computed<DrawerField[]>(() => {
    const stock = props.stock
    if (!stock) return []

    return [
      {
        label: t('inventory.stock_drawer.fields.stock_unit'),
        value: stockUnitLabel.value,
      },
      {
        label: t('inventory.stock_drawer.fields.base_unit'),
        value: baseUnitLabel.value,
      },
      {
        label: t('inventory.stock_drawer.fields.base_units_per_unit'),
        value: toDisplayValue(stock.stock_unit.base_units_per_unit),
      },
    ]
  })

  const materialIdentityFields = computed<DrawerField[]>(() => {
    const stock = props.stock
    if (!stock) return []

    return [
      {
        label: t('inventory.stock_drawer.fields.device_type'),
        value: toDisplayValue(stock.material.device_type?.label || stock.material.device_type?.name),
      },
      {
        label: t('inventory.stock_drawer.fields.item_type'),
        value: toDisplayValue(stock.material.item_type?.label || stock.material.item_type?.name),
      },
      {
        label: t('inventory.stock_drawer.fields.attributes'),
        value: attributeLabel.value,
        wide: true,
      },
      {
        label: t('inventory.stock_drawer.fields.capacity'),
        value: capacityWithUnitLabel.value,
      },
      {
        label: t('inventory.stock_drawer.fields.description'),
        value: toDisplayValue(props.materialDetail?.description),
        wide: true,
      },
    ]
  })

  const supplierFields = computed<DrawerField[]>(() => {
    const stock = props.stock
    if (!stock) return []

    return [
      {
        label: t('inventory.stock_drawer.fields.manufacturer'),
        value: toDisplayValue(stock.material.manufacturer?.label || stock.material.manufacturer?.name),
      },
      {
        label: t('inventory.stock_drawer.fields.vendor'),
        value: toDisplayValue(stock.material.vendor?.label || stock.material.vendor?.name),
      },
      {
        label: t('inventory.stock_drawer.fields.manufacturer_catalog_number'),
        value: toDisplayValue(stock.material.manufacturer_catalog_number),
        wide: true,
      },
      {
        label: t('inventory.stock_drawer.fields.vendor_catalog_number'),
        value: toDisplayValue(stock.material.vendor_catalog_number),
        wide: true,
      },
      {
        label: t('inventory.stock_drawer.fields.brand'),
        value: toDisplayValue(stock.material.brand?.label || stock.material.brand?.name),
      },
      {
        label: t('inventory.stock_drawer.fields.default_cost'),
        value: toDisplayValue(stock.material.default_cost),
      },
    ]
  })

  const sortedUsageEntries = computed<InventoryUsageListItem[]>(() => {
    const usages = [...props.usageEntries]
    usages.sort((leftEntry, rightEntry) => {
      const leftDate = new Date(leftEntry.used_at).getTime()
      const rightDate = new Date(rightEntry.used_at).getTime()
      return rightDate - leftDate
    })
    return usages
  })

  const linkedExperimentProjectLabels = computed<string[]>(() => {
    const labels = new Set<string>()

    for (const usage of props.usageEntries) {
      const projectLabel = usage.project?.label || usage.project?.name || ''
      const experimentLabel = usage.experiment?.label || usage.experiment?.name || ''

      if (projectLabel !== '' && experimentLabel !== '') {
        labels.add(`${projectLabel} / ${experimentLabel}`)
        continue
      }

      if (projectLabel !== '') {
        labels.add(projectLabel)
        continue
      }

      if (experimentLabel !== '') {
        labels.add(experimentLabel)
      }
    }

    return Array.from(labels)
  })

  const sortedOrderEntries = computed<InventoryOrderListItem[]>(() => {
    const orders = [...props.orderEntries]
    orders.sort((leftOrder, rightOrder) => {
      const leftDate = new Date(leftOrder.order_date).getTime()
      const rightDate = new Date(rightOrder.order_date).getTime()
      return rightDate - leftDate
    })
    return orders
  })

  const responsibleUsers = computed<string[]>(() => {
    const labels = new Set<string>()

    for (const order of props.orderEntries) {
      const userLabel = order.who_ordered?.label || order.who_ordered?.full_name || order.who_ordered?.username || ''
      if (userLabel !== '') {
        labels.add(userLabel)
      }
    }

    return Array.from(labels)
  })

  const metadataFields = computed<DrawerField[]>(() => {
    const stock = props.stock
    if (!stock) return []

    return [
      {
        label: t('inventory.stock_drawer.fields.created_at'),
        value: formatDateTime(
          stock.created_at,
          { dateStyle: 'medium', timeStyle: 'short' },
          t('inventory.stock_drawer.values.none'),
        ),
      },
      {
        label: t('inventory.stock_drawer.fields.updated_at'),
        value: formatDateTime(
          stock.updated_at,
          { dateStyle: 'medium', timeStyle: 'short' },
          t('inventory.stock_drawer.values.none'),
        ),
      },
      {
        label: t('inventory.stock_drawer.fields.serial_number'),
        value: toDisplayValue(props.materialDetail?.serial_number),
      },
      {
        label: t('inventory.stock_drawer.fields.order_number'),
        value: toDisplayValue(props.materialDetail?.order_number),
      },
      {
        label: t('inventory.stock_drawer.fields.lifetime_days'),
        value: toDisplayValue(props.materialDetail?.lifetime_days),
      },
    ]
  })

  const toggleMetadata = (): void => {
    isMetadataExpanded.value = !isMetadataExpanded.value
  }

  return {
    inventoryStatusLabel,
    inventoryStatusColor,
    storageTemperatureLabel,
    operationalFields,
    unitConversionFields,
    materialIdentityFields,
    supplierFields,
    sortedUsageEntries,
    linkedExperimentProjectLabels,
    sortedOrderEntries,
    responsibleUsers,
    isMetadataExpanded,
    metadataFields,
    toggleMetadata,
  }
}
