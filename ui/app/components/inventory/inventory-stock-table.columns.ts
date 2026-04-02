import { h } from 'vue'
import type { ColumnDef } from '@tanstack/vue-table'
import type { TableRow } from '~/components/tables/base-data-table.utils'
import type { InventoryStockListItem } from '~/types/inventory'

export type InventoryStockTableRow = TableRow & {
  stock: InventoryStockListItem
  favoriteFlag: boolean
  inventoryStatusCode: InventoryStockListItem['inventory_status']
  inventoryStatusLabel: string
  productName: string
  quantityWithStockUnit: string
  minimumQuantity: string
  location: string
  deviceType: string
  itemType: string
  attributes: string[]
  lotNumber: string
  expiryDate: string
  notes: string
}

type TranslateFn = (key: string, named?: Record<string, unknown>) => string

/**
 * Creates all visible stock table columns used for operational scanning.
 *
 * Returned data example:
 * - `[{ accessorKey: 'favoriteFlag', header: 'Favorite' }, { accessorKey: 'inventoryStatusLabel', header: 'Status' }]`
 */
export const createInventoryStockTableColumns = (t: TranslateFn): ColumnDef<TableRow, unknown>[] => {
  const lowStockBadgeClass =
    'inline-flex items-center rounded-full border border-amber-300 bg-amber-100 px-2 py-0.5 text-[11px] font-semibold text-amber-900'
  const inStockBadgeClass =
    'inline-flex items-center rounded-full border border-emerald-300 bg-emerald-100 px-2 py-0.5 text-[11px] font-semibold text-emerald-900'

  const lifecycleMeta = { semanticGroup: 'lifecycle' } as unknown as Record<string, unknown>
  const identityMeta = { semanticGroup: 'identity' } as unknown as Record<string, unknown>

  const columns: ColumnDef<InventoryStockTableRow, unknown>[] = [
    {
      accessorKey: 'favoriteFlag',
      header: t('inventory.stock_table.columns.favorite'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 78,
      cell: ({ row }) => {
        const isFavorite = row.original.favoriteFlag
        return h(
          'span',
          {
            class: isFavorite ? 'text-amber-500' : 'text-slate-300',
            title: isFavorite ? t('inventory.stock_table.values.yes') : t('inventory.stock_table.values.no'),
          },
          isFavorite ? '★' : '☆',
        )
      },
    },
    {
      accessorKey: 'inventoryStatusLabel',
      header: t('inventory.stock_table.columns.inventory_status'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 132,
      cell: ({ row }) => {
        const isLowStock = row.original.inventoryStatusCode === 'low'
        return h(
          'span',
          {
            class: isLowStock ? lowStockBadgeClass : inStockBadgeClass,
          },
          row.original.inventoryStatusLabel,
        )
      },
    },
    {
      accessorKey: 'productName',
      header: t('inventory.stock_table.columns.product_name'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 260,
    },
    {
      accessorKey: 'quantityWithStockUnit',
      header: t('inventory.stock_table.columns.quantity_with_stock_unit'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 190,
    },
    {
      accessorKey: 'location',
      header: t('inventory.stock_table.columns.location'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 220,
    },
    {
      accessorKey: 'minimumQuantity',
      header: t('inventory.stock_table.columns.minimum_quantity'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 170,
    },
    {
      accessorKey: 'deviceType',
      header: t('inventory.stock_table.columns.device_type'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 180,
      meta: identityMeta,
    },
    {
      accessorKey: 'itemType',
      header: t('inventory.stock_table.columns.item_type'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 180,
      meta: identityMeta,
    },
    {
      accessorKey: 'attributes',
      header: t('inventory.stock_table.columns.attributes'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 220,
      meta: identityMeta,
    },
    {
      accessorKey: 'lotNumber',
      header: t('inventory.stock_table.columns.lot_number'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 170,
      meta: lifecycleMeta,
    },
    {
      accessorKey: 'expiryDate',
      header: t('inventory.stock_table.columns.expiry_date'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 170,
      meta: lifecycleMeta,
    },
    {
      accessorKey: 'notes',
      header: t('inventory.stock_table.columns.notes'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 260,
      meta: lifecycleMeta,
    },
  ]

  return columns as unknown as ColumnDef<TableRow, unknown>[]
}
