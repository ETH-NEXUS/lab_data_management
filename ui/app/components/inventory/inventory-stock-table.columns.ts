import { h } from 'vue'
import type { ColumnDef } from '@tanstack/vue-table'
import type { TableRow } from '~/components/tables/base-data-table.utils'
import type { InventoryStockListItem } from '~/types/inventory'
import {
  getInventoryStockTableSortValue,
  getStatusLabel,
  type TranslateFn,
} from '~/components/inventory/inventory-stock-table.values'

const getStock = (row: TableRow): InventoryStockListItem => {
  return row as InventoryStockListItem
}

/**
 * Creates all visible stock table columns used for operational scanning.
 *
 * Returned data example:
 * - `[{ id: 'productName', header: 'Product' }, { id: 'inventoryStatus', header: 'Status' }]`
 */
export const createInventoryStockTableColumns = (
  t: TranslateFn,
  onSelectStock: (stock: InventoryStockListItem) => void,
  onOpenRelatedOrder: (orderId: number) => void,
): ColumnDef<TableRow, unknown>[] => {
  const outOfStockBadgeClass =
    'inline-flex items-center rounded-full border border-rose-300 bg-rose-100 px-2 py-0.5 text-[11px] font-semibold text-rose-900'
  const lowStockBadgeClass =
    'inline-flex items-center rounded-full border border-amber-300 bg-amber-100 px-2 py-0.5 text-[11px] font-semibold text-amber-900'
  const inStockBadgeClass =
    'inline-flex items-center rounded-full border border-emerald-300 bg-emerald-100 px-2 py-0.5 text-[11px] font-semibold text-emerald-900'

  const lifecycleMeta = { semanticGroup: 'lifecycle' } as unknown as Record<string, unknown>
  const identityMeta = { semanticGroup: 'identity' } as unknown as Record<string, unknown>

  const columns: ColumnDef<TableRow, unknown>[] = [
    {
      id: 'productName',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'productName', t),
      header: t('inventory.stock_table.columns.product_name'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 220,
      cell: ({ row }) => {
        const stock = getStock(row.original)
        const productName = stock.material.product_name?.trim() || t('inventory.stock_table.values.none')
        const label = stock.is_favorite ? `★ ${productName}` : productName
        const linkedOrderId = stock.source_order?.id ?? null

        return h('div', { class: 'flex flex-col items-start gap-1' }, [
          h(
            'button',
            {
              type: 'button',
              class:
                'inline-flex w-full cursor-pointer items-center justify-start rounded-sm px-1 py-0.5 text-left text-blue-700 transition hover:text-blue-800 hover:underline focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-blue-300',
              title: t('inventory.stock_table.columns.product_name'),
              onClick: () => onSelectStock(stock),
            },
            label,
          ),
          linkedOrderId === null
            ? null
            : h(
                'button',
                {
                  type: 'button',
                  class:
                    'inline-flex cursor-pointer items-center rounded-sm px-1 py-0.5 text-xs font-medium text-sky-700 transition hover:text-sky-800 hover:underline focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-blue-300',
                  onClick: () => onOpenRelatedOrder(linkedOrderId),
                },
                `[Order #${linkedOrderId}]`,
              ),
        ])
      },
    },
    // Shows the same favorite flag already visible as a star in the product name,
    // as a plain sortable/filterable "Yes"/"No" column.
    {
      id: 'favorite',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'favorite', t),
      header: t('inventory.stock_table.columns.favorite'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 100,
    },
    {
      id: 'inventoryStatus',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'inventoryStatus', t),
      header: t('inventory.stock_table.columns.inventory_status'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 122,
      cell: ({ row }) => {
        const stock = getStock(row.original)
        const badgeClass =
          stock.inventory_status === 'out_of_stock'
            ? outOfStockBadgeClass
            : stock.inventory_status === 'low'
              ? lowStockBadgeClass
              : inStockBadgeClass
        return h(
          'span',
          {
            class: badgeClass,
          },
          getStatusLabel(t, stock.inventory_status),
        )
      },
    },
    {
      id: 'quantityWithStockUnit',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'quantityWithStockUnit', t),
      header: t('inventory.stock_table.columns.quantity_with_stock_unit'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 150,
    },
    {
      id: 'minimumQuantity',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'minimumQuantity', t),
      header: t('inventory.stock_table.columns.minimum_quantity'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 130,
    },
    {
      id: 'location',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'location', t),
      header: t('inventory.stock_table.columns.location'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 150,
    },
    {
      id: 'deviceType',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'deviceType', t),
      header: t('inventory.stock_table.columns.device_type'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 150,
      meta: identityMeta,
    },
    {
      id: 'itemType',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'itemType', t),
      header: t('inventory.stock_table.columns.item_type'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 150,
      meta: identityMeta,
    },
    // Split out of the itemType cell into its own column so it can be sorted/filtered on its own.
    {
      id: 'storageTemperature',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'storageTemperature', t),
      header: t('inventory.stock_table.columns.storage_temperature'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 140,
      meta: identityMeta,
    },
    {
      id: 'attributes',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'attributes', t),
      header: t('inventory.stock_table.columns.attributes'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 180,
      meta: identityMeta,
    },
    // Batch 1: material identity fields already shown in the stock detail drawer's
    // "Supplier and catalog" section, now also available directly in the table.
    {
      id: 'brand',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'brand', t),
      header: t('inventory.stock_table.columns.brand'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 140,
      meta: identityMeta,
    },
    {
      id: 'manufacturer',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'manufacturer', t),
      header: t('inventory.stock_table.columns.manufacturer'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 150,
      meta: identityMeta,
    },
    {
      id: 'vendor',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'vendor', t),
      header: t('inventory.stock_table.columns.vendor'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 150,
      meta: identityMeta,
    },
    // Batch 2: catalog references and remaining supplier/catalog fields from the
    // stock detail drawer's "Supplier and catalog" section.
    {
      id: 'manufacturerCatalogNumber',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'manufacturerCatalogNumber', t),
      header: t('inventory.stock_table.columns.manufacturer_catalog_number'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 180,
      meta: identityMeta,
    },
    {
      id: 'vendorCatalogNumber',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'vendorCatalogNumber', t),
      header: t('inventory.stock_table.columns.vendor_catalog_number'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 180,
      meta: identityMeta,
    },
    {
      id: 'capacity',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'capacity', t),
      header: t('inventory.stock_table.columns.capacity'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 140,
      meta: identityMeta,
    },
    {
      id: 'defaultCost',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'defaultCost', t),
      header: t('inventory.stock_table.columns.default_cost'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 130,
      meta: identityMeta,
    },
    {
      id: 'isActive',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'isActive', t),
      header: t('inventory.stock_table.columns.is_active'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 110,
      meta: identityMeta,
    },
    {
      id: 'lotNumber',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'lotNumber', t),
      header: t('inventory.stock_table.columns.lot_number'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 130,
      meta: lifecycleMeta,
    },
    {
      id: 'expiryDate',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'expiryDate', t),
      header: t('inventory.stock_table.columns.expiry_date'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 130,
      meta: lifecycleMeta,
    },
    {
      id: 'notes',
      accessorFn: (row) => getInventoryStockTableSortValue(getStock(row), 'notes', t),
      header: t('inventory.stock_table.columns.notes'),
      enableSorting: true,
      enableColumnFilter: true,
      size: 210,
      meta: lifecycleMeta,
    },
  ]

  return columns
}
