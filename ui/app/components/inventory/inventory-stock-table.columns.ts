import { h } from 'vue'
import type { ColumnDef } from '@tanstack/vue-table'
import type { TableRow } from '~/components/tables/base-data-table.utils'
import type { InventoryStockListItem } from '~/types/inventory'
import { formatDateTime } from '~/utils/dateTime'

type TranslateFn = (key: string, named?: Record<string, unknown>) => string

const toLabel = (value: string | null | undefined, fallback: string): string => {
  return value?.trim() || fallback
}

const getStock = (row: TableRow): InventoryStockListItem => {
  return row as InventoryStockListItem
}

const getStatusLabel = (t: TranslateFn, status: InventoryStockListItem['inventory_status']): string => {
  return status === 'low'
    ? t('inventory.stock_table.status_labels.low')
    : t('inventory.stock_table.status_labels.in_stock')
}

const getAttributes = (stock: InventoryStockListItem, t: TranslateFn): string[] => {
  const labels = (stock.material.attributes ?? [])
    .map((attribute) => attribute.label || attribute.name || '')
    .map((label) => label.trim())
    .filter((label) => label !== '')

  return labels.length > 0 ? labels : [t('inventory.stock_drawer.values.no_attributes')]
}

/**
 * Creates all visible stock table columns used for operational scanning.
 *
 * Returned data example:
 * - `[{ id: 'productName', header: 'Product' }, { id: 'inventoryStatus', header: 'Status' }]`
 */
export const createInventoryStockTableColumns = (t: TranslateFn): ColumnDef<TableRow, unknown>[] => {
  const lowStockBadgeClass =
    'inline-flex items-center rounded-full border border-amber-300 bg-amber-100 px-2 py-0.5 text-[11px] font-semibold text-amber-900'
  const inStockBadgeClass =
    'inline-flex items-center rounded-full border border-emerald-300 bg-emerald-100 px-2 py-0.5 text-[11px] font-semibold text-emerald-900'

  const lifecycleMeta = { semanticGroup: 'lifecycle' } as unknown as Record<string, unknown>
  const identityMeta = { semanticGroup: 'identity' } as unknown as Record<string, unknown>

  const columns: ColumnDef<TableRow, unknown>[] = [
    {
      id: 'productName',
      accessorFn: (row) => toLabel(getStock(row).material.product_name, t('inventory.stock_table.values.none')),
      header: t('inventory.stock_table.columns.product_name'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 260,
      cell: ({ row }) => {
        const stock = getStock(row.original)
        const productName = toLabel(stock.material.product_name, t('inventory.stock_table.values.none'))
        if (!stock.is_favorite) {
          return productName
        }

        return h(
          'span',
          {
            class:
              'block w-full rounded-md border border-amber-300 bg-amber-100 px-2 py-1 font-semibold text-amber-900',
            title: t('inventory.stock_table.values.yes'),
          },
          `★ ${productName}`,
        )
      },
    },
    {
      id: 'inventoryStatus',
      accessorFn: (row) => getStatusLabel(t, getStock(row).inventory_status),
      header: t('inventory.stock_table.columns.inventory_status'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 132,
      cell: ({ row }) => {
        const stock = getStock(row.original)
        const isLowStock = stock.inventory_status === 'low'
        return h(
          'span',
          {
            class: isLowStock ? lowStockBadgeClass : inStockBadgeClass,
          },
          getStatusLabel(t, stock.inventory_status),
        )
      },
    },
    {
      id: 'quantityWithStockUnit',
      accessorFn: (row) => toLabel(getStock(row).stock_label, t('inventory.stock_table.values.none')),
      header: t('inventory.stock_table.columns.quantity_with_stock_unit'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 190,
    },
    {
      id: 'minimumQuantity',
      accessorFn: (row) => toLabel(getStock(row).minimum_quantity, t('inventory.stock_table.values.none')),
      header: t('inventory.stock_table.columns.minimum_quantity'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 170,
    },
    {
      id: 'location',
      accessorFn: (row) => toLabel(getStock(row).location_label, t('inventory.stock_table.values.unknown_location')),
      header: t('inventory.stock_table.columns.location'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 220,
    },
    {
      id: 'deviceType',
      accessorFn: (row) => {
        const stock = getStock(row)
        return toLabel(
          stock.material.device_type?.label || stock.material.device_type?.name,
          t('inventory.stock_table.values.none'),
        )
      },
      header: t('inventory.stock_table.columns.device_type'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 180,
      meta: identityMeta,
    },
    {
      id: 'itemType',
      accessorFn: (row) => {
        const stock = getStock(row)
        return toLabel(
          stock.material.item_type?.label || stock.material.item_type?.name,
          t('inventory.stock_table.values.none'),
        )
      },
      header: t('inventory.stock_table.columns.item_type'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 180,
      meta: identityMeta,
    },
    {
      id: 'attributes',
      accessorFn: (row) => getAttributes(getStock(row), t),
      header: t('inventory.stock_table.columns.attributes'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 220,
      meta: identityMeta,
    },
    {
      id: 'lotNumber',
      accessorFn: (row) => toLabel(getStock(row).lot_number, t('inventory.stock_table.values.none')),
      header: t('inventory.stock_table.columns.lot_number'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 170,
      meta: lifecycleMeta,
    },
    {
      id: 'expiryDate',
      accessorFn: (row) => {
        const expiryDate = getStock(row).expiry_date
        if (!expiryDate) {
          return t('inventory.stock_table.values.none')
        }
        return formatDateTime(expiryDate, { dateStyle: 'medium' }, t('inventory.stock_table.values.none'))
      },
      header: t('inventory.stock_table.columns.expiry_date'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 170,
      meta: lifecycleMeta,
    },
    {
      id: 'notes',
      accessorFn: (row) => toLabel(getStock(row).notes, t('inventory.stock_table.values.no_notes')),
      header: t('inventory.stock_table.columns.notes'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 260,
      meta: lifecycleMeta,
    },
  ]

  return columns
}
