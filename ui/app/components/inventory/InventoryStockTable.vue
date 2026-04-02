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
 * Resolves location label using API-provided combined value.
 *
 * Accepted stock example:
 * - `{ location_label: 'Room A / Shelf 3' }`
 *
 * Returned data examples:
 * - `'Room A / Shelf 3'`
 * - `'Unknown location'`
 */
const formatLocationLabel = (stock: InventoryStockListItem): string => {
  return stock.location_label?.trim() || t('inventory.stock_table.values.unknown_location')
}

/**
 * Builds quantity label using API-provided stock label.
 *
 * Input example:
 * - `{ stock_label: '2.5 L' }`
 *
 * Returned example:
 * - `'2.5 L'`
 */
const formatQuantityWithStockUnit = (stock: InventoryStockListItem): string => {
  return stock.stock_label?.trim() || t('inventory.stock_table.values.none')
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
  return props.stocks.map((stock) => {
    const attributeLabels = (stock.material.attributes ?? [])
      .map((attribute) => attribute.label || attribute.name || '')
      .map((label) => label.trim())
      .filter((label) => label !== '')

    return {
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
    }
  })
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
