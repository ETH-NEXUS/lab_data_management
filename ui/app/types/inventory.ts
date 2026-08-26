export const INVENTORY_STOCKS_ENDPOINT = 'inventory/stocks/'
export const INVENTORY_MATERIALS_ENDPOINT = 'inventory/materials/'
export const INVENTORY_MATERIAL_UNITS_ENDPOINT = 'inventory/material-units/'
export const INVENTORY_ORDERS_ENDPOINT = 'inventory/orders/'
export const INVENTORY_USAGES_ENDPOINT = 'inventory/material-usages/'
export const INVENTORY_STOCK_TABLE_PREFERENCES_ENDPOINT = 'inventory/stock-table-preferences/'

export const INVENTORY_ROOMS_ENDPOINT = 'inventory/rooms/'
export const INVENTORY_SECTORS_ENDPOINT = 'inventory/sectors/'
export const INVENTORY_MANUFACTURERS_ENDPOINT = 'inventory/manufacturers/'
export const INVENTORY_BRANDS_ENDPOINT = 'inventory/brands/'
export const INVENTORY_VENDORS_ENDPOINT = 'inventory/vendors/'
export const INVENTORY_DEVICE_TYPES_ENDPOINT = 'inventory/device-types/'
export const INVENTORY_ITEM_TYPES_ENDPOINT = 'inventory/item-types/'
export const INVENTORY_ATTRIBUTES_ENDPOINT = 'inventory/material-attributes/'
export const INVENTORY_UNITS_ENDPOINT = 'inventory/units-of-measure/'

export const INVENTORY_STOCKS_QUERY_KEY = ['inventory-stocks']
export const INVENTORY_STOCK_PAGES_QUERY_KEY = ['inventory-stocks-page']
export const INVENTORY_MATERIALS_QUERY_KEY = ['inventory-materials']
export const INVENTORY_MATERIAL_UNITS_QUERY_KEY = ['inventory-material-units']
export const INVENTORY_ORDERS_QUERY_KEY = ['inventory-orders']
export const INVENTORY_USAGES_QUERY_KEY = ['inventory-usages']
export const INVENTORY_LOOKUPS_QUERY_KEY = ['inventory-lookups']
export const INVENTORY_STOCK_TABLE_PREFERENCES_QUERY_KEY = ['inventory-stock-table-preferences']

export const INVENTORY_STOCKS_ERROR_MESSAGE = 'Failed to load inventory stock items.'
export const INVENTORY_STOCK_ERROR_MESSAGE = 'Failed to load inventory stock item.'
export const INVENTORY_MATERIALS_ERROR_MESSAGE = 'Failed to load inventory materials.'
export const INVENTORY_MATERIAL_ERROR_MESSAGE = 'Failed to load inventory material.'
export const INVENTORY_MATERIAL_UNITS_ERROR_MESSAGE = 'Failed to load inventory material units.'
export const INVENTORY_ORDERS_ERROR_MESSAGE = 'Failed to load inventory orders.'
export const INVENTORY_ORDER_ERROR_MESSAGE = 'Failed to load inventory order.'
export const INVENTORY_USAGES_ERROR_MESSAGE = 'Failed to load material usages.'
export const INVENTORY_USAGE_ERROR_MESSAGE = 'Failed to load material usage.'
export const INVENTORY_CREATE_STOCK_ERROR_MESSAGE = 'Failed to create inventory stock item.'
export const INVENTORY_UPDATE_STOCK_ERROR_MESSAGE = 'Failed to update inventory stock item.'
export const INVENTORY_CREATE_MATERIAL_ERROR_MESSAGE = 'Failed to create inventory material.'
export const INVENTORY_UPDATE_MATERIAL_ERROR_MESSAGE = 'Failed to update inventory material.'
export const INVENTORY_CREATE_ORDER_ERROR_MESSAGE = 'Failed to create order.'
export const INVENTORY_UPDATE_ORDER_ERROR_MESSAGE = 'Failed to update order.'
export const INVENTORY_CREATE_USAGE_ERROR_MESSAGE = 'Failed to create material usage.'
export const INVENTORY_UPDATE_USAGE_ERROR_MESSAGE = 'Failed to update material usage.'
export const INVENTORY_MARK_FAVORITE_ERROR_MESSAGE = 'Failed to mark stock item as favorite.'
export const INVENTORY_UNMARK_FAVORITE_ERROR_MESSAGE = 'Failed to unmark stock item as favorite.'
export const INVENTORY_ARCHIVE_STOCK_ERROR_MESSAGE = 'Failed to archive stock item.'
export const INVENTORY_RESTORE_STOCK_ERROR_MESSAGE = 'Failed to restore stock item.'
export const INVENTORY_STOCK_TABLE_PREFERENCE_ERROR_MESSAGE = 'Failed to load inventory stock table preference.'
export const INVENTORY_UPDATE_STOCK_TABLE_PREFERENCE_ERROR_MESSAGE =
  'Failed to update inventory stock table preference.'

export const getInventoryStockQueryKey = (stockId: number) => [...INVENTORY_STOCKS_QUERY_KEY, stockId]
export const getInventoryMaterialQueryKey = (materialId: number) => [...INVENTORY_MATERIALS_QUERY_KEY, materialId]
export const getInventoryOrderQueryKey = (orderId: number) => [...INVENTORY_ORDERS_QUERY_KEY, orderId]
export const getInventoryUsageQueryKey = (usageId: number) => [...INVENTORY_USAGES_QUERY_KEY, usageId]
export const getInventoryStockTablePreferenceQueryKey = (tableKey: string) => [
  ...INVENTORY_STOCK_TABLE_PREFERENCES_QUERY_KEY,
  tableKey,
]

export type InventoryLookupItem = {
  id: number
  name: string
  label: string
}

export type InventoryRoom = {
  id: number
  name: string
  label: string
}

export type InventorySector = {
  id: number
  name: string
  label: string
  room: InventoryRoom
}

export type InventoryMaterialAttribute = {
  id: number
  name: string
  label: string
}

export type InventoryUnitOfMeasure = {
  id: number
  name: string
  label: string
}

export type InventoryMaterialUnitSummary = {
  id: number
  unit: InventoryUnitOfMeasure
  is_base_unit: boolean
  is_stock_unit: boolean
  is_order_unit: boolean
  base_units_per_unit: string
  display_name: string
}

export type InventoryMaterialListItem = {
  id: number
  product_name: string
  label: string
  brand: InventoryLookupItem | null
  manufacturer: InventoryLookupItem | null
  vendor: InventoryLookupItem | null
  manufacturer_catalog_number: string | null
  vendor_catalog_number: string | null
  item_type: InventoryLookupItem | null
  device_type: InventoryLookupItem | null
  attributes: InventoryMaterialAttribute[]
  capacity_value: string | null
  capacity_unit: string | null
  capacity_display: string | null
  storage_temperature: string | null
  storage_temperature_label: string | null
  safety_data_sheet: string | null
  default_cost: string | null
  description: string | null
  lifetime_days: number | null
  serial_number: string | null
  order_number: string | null
  is_active: boolean
}

export type InventoryMaterialDetail = InventoryMaterialListItem & {
  units: InventoryMaterialUnitSummary[]
}

export type InventoryOrderSourceSummary = {
  id: number
  order_date: string
  status: string
  status_label: string
  created_at: string
}

export type InventoryCreatedStockSummary = {
  id: number
  lot_number: string | null
  expiry_date: string | null
  location_label: string | null
  created_at: string
}

export type InventoryStockListItem = {
  id: number
  material: InventoryMaterialListItem
  sector: InventorySector
  sectors: InventorySector[]
  room: InventoryRoom | null
  stock_unit: InventoryMaterialUnitSummary
  quantity: string
  minimum_quantity: string
  inventory_status: 'out_of_stock' | 'low' | 'in stock'
  quantity_in_base_units: string
  minimum_quantity_in_base_units: string
  lot_number: string | null
  expiry_date: string | null
  is_favorite: boolean
  is_archived: boolean
  archived_at: string | null
  notes: string | null
  source_order: InventoryOrderSourceSummary | null
  location_label: string | null
  stock_label: string | null
  is_low_stock: boolean
  is_expired: boolean
  is_expiring_soon: boolean
  created_at: string
  updated_at: string
}

export type InventoryStockDetail = InventoryStockListItem & {
  material_id?: number
  sector_id?: number
  sector_ids?: number[]
  stock_unit_id?: number
}

export type InventoryUserSummary = {
  id: number
  username: string
  full_name: string
  label: string
}

export type InventoryProjectSummary = {
  id: number
  name?: string
  label?: string
}

export type InventoryExperimentSummary = {
  id: number
  name?: string
  label?: string
}

export type InventoryOrderListItem = {
  id: number
  material: InventoryMaterialListItem
  order_unit: InventoryMaterialUnitSummary
  amount: string
  order_date: string
  status: string
  status_label: string
  who_ordered: InventoryUserSummary | null
  project: InventoryProjectSummary | null
  notes: string | null
  created_stock_entries: InventoryCreatedStockSummary[]
  created_at: string
  updated_at: string
}

export type InventoryOrderDetail = InventoryOrderListItem & {
  material_id?: number
  order_unit_id?: number
  who_ordered_id?: number | null
  project_id?: number | null
}

export type InventoryUsageListItem = {
  id: number
  project: InventoryProjectSummary
  experiment: InventoryExperimentSummary | null
  inventory_stock: InventoryStockListItem
  material: InventoryMaterialListItem | null
  usage_unit: InventoryMaterialUnitSummary
  quantity_used: string
  quantity_used_in_base_units: string
  used_at: string
  notes: string | null
}

export type InventoryUsageDetail = InventoryUsageListItem & {
  project_id?: number
  experiment_id?: number | null
  inventory_stock_id?: number
  usage_unit_id?: number
}

export type InventoryStockPreset = 'all' | 'favorite' | 'low_stock' | 'expired' | 'archived'

export type InventoryStockTableSorting = {
  id: string
  desc: boolean
}

export type InventoryStockTableColumnFilter = {
  id: string
  value: string[]
}

export type InventoryStockTablePreference = {
  id: number
  table_key: string
  preset: InventoryStockPreset
  sorting: InventoryStockTableSorting[]
  column_filters: InventoryStockTableColumnFilter[]
  column_order: string[]
  column_visibility: Record<string, boolean>
  created_at: string
  updated_at: string
}

export type CreateInventoryMaterialPayload = {
  product_name: string
  brand_id?: number | null
  manufacturer_id?: number | null
  vendor_id?: number | null
  manufacturer_catalog_number?: string | null
  vendor_catalog_number?: string | null
  item_type_id?: number | null
  device_type_id?: number | null
  attribute_ids?: number[]
  capacity_value?: string | null
  capacity_unit?: string | null
  storage_temperature?: string | null
  safety_data_sheet?: File | null
  description?: string | null
  default_cost?: string | null
  lifetime_days?: number | null
  serial_number?: string | null
  order_number?: string | null
  is_active?: boolean
}

export type UpdateInventoryMaterialPayload = Partial<CreateInventoryMaterialPayload>

export type CreateInventoryStockPayload = {
  material_id: number
  sector_ids: number[]
  stock_unit_id: number
  quantity: string
  minimum_quantity: string
  source_order_id?: number | null
  lot_number?: string | null
  expiry_date?: string | null
  notes?: string | null
}

export type UpdateInventoryStockPayload = Partial<CreateInventoryStockPayload>

export type CreateInventoryOrderPayload = {
  material_id: number
  order_unit_id: number
  amount: string
  order_date: string
  status: string
  who_ordered_id?: number | null
  project_id?: number | null
  notes?: string | null
}

export type UpdateInventoryOrderPayload = Partial<CreateInventoryOrderPayload>

export type CreateInventoryUsagePayload = {
  project_id: number
  experiment_id?: number | null
  inventory_stock_id: number
  usage_unit_id: number
  quantity_used: string
  notes?: string | null
}

export type UpdateInventoryUsagePayload = Partial<CreateInventoryUsagePayload>
