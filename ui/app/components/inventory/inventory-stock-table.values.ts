import type { InventoryStockListItem, InventoryStockPreset } from '~/types/inventory'

export type TranslateFn = (key: string, named?: Record<string, unknown>) => string

const toLabel = (value: string | null | undefined, fallback: string): string => {
  return value?.trim() || fallback
}

const formatNumericString = (value: string | null | undefined): string => {
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

  return rawValue.replace(/\.?0+$/, '')
}

export const getStatusLabel = (t: TranslateFn, status: InventoryStockListItem['inventory_status']): string => {
  if (status === 'out_of_stock') {
    return t('inventory.stock_table.status_labels.out_of_stock')
  }

  return status === 'low'
    ? t('inventory.stock_table.status_labels.low')
    : t('inventory.stock_table.status_labels.in_stock')
}

export const getAttributes = (stock: InventoryStockListItem, t: TranslateFn): string[] => {
  const labels = (stock.material.attributes ?? [])
    .map((attribute) => attribute.label || attribute.name || '')
    .map((label) => label.trim())
    .filter((label) => label !== '')

  return labels.length > 0 ? labels : [t('inventory.stock_drawer.values.no_attributes')]
}

export const getStockUnitLabel = (stock: InventoryStockListItem): string => {
  return stock.stock_unit.unit?.label || stock.stock_unit.unit?.name || ''
}

export const getFormattedQuantityWithUnit = (stock: InventoryStockListItem, fallback: string): string => {
  const quantityValue = formatNumericString(stock.quantity)
  if (quantityValue === '') {
    return fallback
  }

  const unitLabel = getStockUnitLabel(stock)
  if (unitLabel === '') {
    return quantityValue
  }

  return `${quantityValue} ${unitLabel}`.trim()
}

/**
 * Applies one dashboard/workspace preset to an already-loaded stock list.
 *
 * Input examples:
 * - `preset = 'favorite'`
 * - `preset = 'low_stock'`
 *
 * Returned data example:
 * - `[{ id: 4, is_favorite: true, ... }, { id: 12, is_favorite: true, ... }]`
 */
export const getStocksForPreset = (
  stocks: InventoryStockListItem[],
  preset: InventoryStockPreset,
): InventoryStockListItem[] => {
  if (preset === 'favorite') {
    return stocks.filter((stock) => stock.is_favorite)
  }

  if (preset === 'low_stock') {
    return stocks.filter((stock) => stock.is_low_stock || stock.inventory_status === 'low')
  }

  if (preset === 'expired') {
    return stocks.filter((stock) => stock.is_expired)
  }

  return stocks
}

/**
 * Returns the same accessor-style value shape used by the stock table for one sortable column.
 *
 * Input examples:
 * - `columnId = 'productName'`
 * - `columnId = 'expiryDate'`
 *
 * Returned data examples:
 * - `'PBS Buffer'`
 * - `'Low stock'`
 * - `['filtered', 'sterile']`
 */
export const getInventoryStockTableSortValue = (
  stock: InventoryStockListItem,
  columnId: string,
  t: TranslateFn,
): string | string[] => {
  if (columnId === 'productName') {
    return toLabel(stock.material.product_name, t('inventory.stock_table.values.none'))
  }

  // Accepted data example: `{ is_favorite: true }`
  // Returned data example: `'Yes'`
  if (columnId === 'favorite') {
    return stock.is_favorite ? t('inventory.stock_table.values.yes') : t('inventory.stock_table.values.no')
  }

  if (columnId === 'inventoryStatus') {
    return getStatusLabel(t, stock.inventory_status)
  }

  if (columnId === 'quantityWithStockUnit') {
    return getFormattedQuantityWithUnit(stock, t('inventory.stock_table.values.none'))
  }

  if (columnId === 'minimumQuantity') {
    return toLabel(formatNumericString(stock.minimum_quantity), t('inventory.stock_table.values.none'))
  }

  if (columnId === 'location') {
    return toLabel(stock.location_label, t('inventory.stock_table.values.unknown_location'))
  }

  if (columnId === 'deviceType') {
    return toLabel(
      stock.material.device_type?.label || stock.material.device_type?.name,
      t('inventory.stock_table.values.none'),
    )
  }

  if (columnId === 'itemType') {
    return toLabel(
      stock.material.item_type?.label || stock.material.item_type?.name,
      t('inventory.stock_table.values.none'),
    )
  }

  // Accepted data example: `{ storage_temperature_label: '4°C', storage_temperature: '4C' }`
  // Returned data example: `'4°C'`
  if (columnId === 'storageTemperature') {
    return toLabel(
      stock.material.storage_temperature_label || stock.material.storage_temperature,
      t('inventory.stock_table.values.none'),
    )
  }

  if (columnId === 'attributes') {
    return getAttributes(stock, t)
  }

  // Accepted data example: `{ brand: { name: 'Costar', label: 'Costar' } }`
  // Returned data example: `'Costar'`
  if (columnId === 'brand') {
    return toLabel(stock.material.brand?.label || stock.material.brand?.name, t('inventory.stock_table.values.none'))
  }

  // Accepted data example: `{ manufacturer: { name: 'Corning', label: 'Corning' } }`
  // Returned data example: `'Corning'`
  if (columnId === 'manufacturer') {
    return toLabel(
      stock.material.manufacturer?.label || stock.material.manufacturer?.name,
      t('inventory.stock_table.values.none'),
    )
  }

  // Accepted data example: `{ vendor: { name: 'Huberlab', label: 'Huberlab' } }`
  // Returned data example: `'Huberlab'`
  if (columnId === 'vendor') {
    return toLabel(
      stock.material.vendor?.label || stock.material.vendor?.name,
      t('inventory.stock_table.values.none'),
    )
  }

  // Accepted data example: `{ manufacturer_catalog_number: '30038616' }`
  // Returned data example: `'30038616'`
  if (columnId === 'manufacturerCatalogNumber') {
    return toLabel(stock.material.manufacturer_catalog_number, t('inventory.stock_table.values.none'))
  }

  // Accepted data example: `{ vendor_catalog_number: '3815' }`
  // Returned data example: `'3815'`
  if (columnId === 'vendorCatalogNumber') {
    return toLabel(stock.material.vendor_catalog_number, t('inventory.stock_table.values.none'))
  }

  // Accepted data example: `{ capacity_display: '200 ul' }`
  // Returned data example: `'200 ul'`
  if (columnId === 'capacity') {
    return toLabel(stock.material.capacity_display, t('inventory.stock_table.values.none'))
  }

  // Accepted data example: `{ default_cost: '12.50' }`
  // Returned data example: `'12.5'`
  if (columnId === 'defaultCost') {
    return toLabel(formatNumericString(stock.material.default_cost), t('inventory.stock_table.values.none'))
  }

  // Accepted data example: `{ is_active: true }`
  // Returned data example: `'Yes'`
  if (columnId === 'isActive') {
    return stock.material.is_active ? t('inventory.stock_table.values.yes') : t('inventory.stock_table.values.no')
  }

  // Accepted data example: `{ description: 'Sterile PCR-clean tips' }`
  // Returned data example: `'Sterile PCR-clean tips'`
  if (columnId === 'description') {
    return toLabel(stock.material.description, t('inventory.stock_table.values.none'))
  }

  // Accepted data example: `{ serial_number: 'SN-4821' }`
  // Returned data example: `'SN-4821'`
  if (columnId === 'serialNumber') {
    return toLabel(stock.material.serial_number, t('inventory.stock_table.values.none'))
  }

  // Accepted data example: `{ order_number: 'PO-2026-0113' }`
  // Returned data example: `'PO-2026-0113'`
  if (columnId === 'orderNumber') {
    return toLabel(stock.material.order_number, t('inventory.stock_table.values.none'))
  }

  // Accepted data example: `{ lifetime_days: 365 }`
  // Returned data example: `'365'`
  if (columnId === 'lifetimeDays') {
    return stock.material.lifetime_days == null
      ? t('inventory.stock_table.values.none')
      : String(stock.material.lifetime_days)
  }

  if (columnId === 'lotNumber') {
    return toLabel(stock.lot_number, t('inventory.stock_table.values.none'))
  }

  if (columnId === 'expiryDate') {
    if (!stock.expiry_date) {
      return t('inventory.stock_table.values.none')
    }

    return stock.expiry_date
  }

  if (columnId === 'notes') {
    return toLabel(stock.notes, t('inventory.stock_table.values.no_notes'))
  }

  return String(stock.id)
}

const normalizeSortValue = (value: string | string[]): string => {
  return Array.isArray(value) ? value.join(', ') : value
}

/**
 * Applies the saved stock table sorting state to one stock list.
 *
 * Input examples:
 * - `[{ id: 'productName', desc: false }]`
 * - `[{ id: 'inventoryStatus', desc: true }, { id: 'productName', desc: false }]`
 */
export const sortStocksLikeInventoryTable = (
  stocks: InventoryStockListItem[],
  sortingState: { id: string; desc: boolean }[],
  t: TranslateFn,
): InventoryStockListItem[] => {
  if (sortingState.length === 0) {
    return stocks
  }

  const sortableStocks = [...stocks]

  sortableStocks.sort((leftStock, rightStock) => {
    for (const sortRule of sortingState) {
      const leftValue = normalizeSortValue(getInventoryStockTableSortValue(leftStock, sortRule.id, t))
      const rightValue = normalizeSortValue(getInventoryStockTableSortValue(rightStock, sortRule.id, t))
      const comparison = leftValue.localeCompare(rightValue, undefined, { numeric: false, sensitivity: 'base' })

      if (comparison !== 0) {
        return sortRule.desc ? comparison * -1 : comparison
      }
    }

    return 0
  })

  return sortableStocks
}
