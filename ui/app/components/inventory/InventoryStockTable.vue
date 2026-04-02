<script setup lang="ts">
import { computed } from 'vue'
import BaseDataTable from '~/components/tables/BaseDataTable.vue'
import type { TableRow } from '~/components/tables/base-data-table.utils'
import {
  createInventoryStockTableColumns,
  type InventoryStockTableRow,
} from '~/components/inventory/inventory-stock-table.columns'
import type { InventoryStockListItem } from '~/types/inventory'
import { formatDateTime } from '~/utils/dateTime'

type Props = {
  stocks: InventoryStockListItem[]
}

const props = defineProps<Props>()

const emit = defineEmits<{
  (e: 'select-stock', stock: InventoryStockListItem): void
}>()

const { t } = useI18n()

/**
 * Resolves one location label from room + sector values.
 *
 * Accepted stock example:
 * - `{ room: { label: 'Room A' }, sector: { label: 'Shelf 3' } }`
 *
 * Returned data examples:
 * - `'Room A / Shelf 3'`
 * - `'Unknown location'`
 */
const formatLocationLabel = (stock: InventoryStockListItem): string => {
  const roomLabel = stock.room?.label ?? stock.room?.name ?? stock.sector?.room?.label ?? stock.sector?.room?.name ?? ''
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

  return t('inventory.stock_table.values.unknown_location')
}

/**
 * Builds one quantity label that includes stock unit name.
 *
 * Input example:
 * - `{ quantity: '2.5', stock_unit: { display_name: 'L' } }`
 *
 * Returned example:
 * - `'2.5 L'`
 */
const formatQuantityWithStockUnit = (stock: InventoryStockListItem): string => {
  const quantityValue = stock.quantity?.trim() || t('inventory.stock_table.values.none')
  const stockUnitLabel =
    stock.stock_unit?.display_name || stock.stock_unit?.unit?.label || stock.stock_unit?.unit?.name || ''

  if (stockUnitLabel === '') {
    return quantityValue
  }

  return `${quantityValue} ${stockUnitLabel}`
}

/**
 * Returns user-visible label for inventory status values.
 *
 * Input examples:
 * - `'in stock'`
 * - `'low'`
 *
 * Returned examples:
 * - `'In stock'`
 * - `'Low stock'`
 */
const getInventoryStatusLabel = (status: InventoryStockListItem['inventory_status']): string => {
  if (status === 'low') {
    return t('inventory.stock_table.status_labels.low')
  }
  return t('inventory.stock_table.status_labels.in_stock')
}

/**
 * Converts raw API stock items to flat table rows used by BaseDataTable.
 *
 * Accepted data example:
 * - `[{ id: 1, material: { product_name: 'PBS' }, quantity: '3', ... }]`
 *
 * Returned data example:
 * - `[{ productName: 'PBS', quantityWithStockUnit: '3 mL', ... }]`
 */
const tableRows = computed<InventoryStockTableRow[]>(() => {
  const rows: InventoryStockTableRow[] = []

  for (const stock of props.stocks) {
    const attributeLabels: string[] = []

    for (const attribute of stock.material.attributes ?? []) {
      const label = attribute.label || attribute.name
      if (label && label.trim() !== '') {
        attributeLabels.push(label)
      }
    }

    rows.push({
      stock,
      favoriteFlag: stock.is_favorite,
      inventoryStatusCode: stock.inventory_status,
      inventoryStatusLabel: getInventoryStatusLabel(stock.inventory_status),
      productName: stock.material.product_name || t('inventory.stock_table.values.none'),
      quantityWithStockUnit: formatQuantityWithStockUnit(stock),
      minimumQuantity: stock.minimum_quantity || t('inventory.stock_table.values.none'),
      location: formatLocationLabel(stock),
      deviceType:
        stock.material.device_type?.label || stock.material.device_type?.name || t('inventory.stock_table.values.none'),
      itemType:
        stock.material.item_type?.label || stock.material.item_type?.name || t('inventory.stock_table.values.none'),
      attributes: attributeLabels.length > 0 ? attributeLabels : [t('inventory.stock_drawer.values.no_attributes')],
      lotNumber: stock.lot_number || t('inventory.stock_table.values.none'),
      expiryDate: stock.expiry_date
        ? formatDateTime(stock.expiry_date, { dateStyle: 'medium' }, t('inventory.stock_table.values.none'))
        : t('inventory.stock_table.values.none'),
      notes: stock.notes || t('inventory.stock_table.values.no_notes'),
    })
  }

  return rows
})

const tableColumns = computed(() => createInventoryStockTableColumns(t))

/**
 * Forwards one selected table row as its original stock payload.
 *
 * Input example:
 * - `row = { stock: { id: 14, ... } }`
 */
const onRowClick = (row: TableRow): void => {
  const tableRow = row as InventoryStockTableRow
  emit('select-stock', tableRow.stock)
}
</script>

<template>
  <BaseDataTable
    :data="tableRows as unknown as TableRow[]"
    :columns="tableColumns"
    :frozen-column-count="5"
    dense
    row-clickable
    :global-filter-placeholder="t('table.general.search_placeholder')"
    @row-click="onRowClick"
  />
</template>
