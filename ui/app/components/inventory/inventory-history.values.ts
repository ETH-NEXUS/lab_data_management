import { formatNumericString } from '~/components/inventory/inventory-stock-table.values'
import type { InventoryHistoryListItem } from '~/types/inventory'

const additionActions = new Set(['stock_created', 'order_created', 'stock_restored'])
const deletionActions = new Set(['stock_deleted', 'order_deleted', 'usage_deleted', 'stock_archived'])

/**
 * Returns a muted row color for the type of change.
 *
 * Examples:
 * - `stock_created` returns a green background.
 * - `usage_created` returns an amber background.
 * - `stock_archived` returns a red background.
 */
export const getHistoryRecordRowClass = (record: InventoryHistoryListItem): string => {
  if (record.performed_action === 'usage_created') {
    return 'bg-[#F6E7C8] hover:bg-[#F1DDB8]'
  }

  if (additionActions.has(record.performed_action)) {
    return 'bg-[#DDEBD8] hover:bg-[#D1E3CB]'
  }

  if (deletionActions.has(record.performed_action)) {
    return 'bg-[#F5DBD5] hover:bg-[#EFCBC3]'
  }

  return 'bg-[#DCEAF7] hover:bg-[#D1E3F2]'
}

export const getHistoryRecordActionLabel = (record: InventoryHistoryListItem): string => {
  return record.performed_action
    .split('_')
    .map((word) => `${word.charAt(0).toUpperCase()}${word.slice(1)}`)
    .join(' ')
}

export const getHistoryRecordItemName = (record: InventoryHistoryListItem, fallback: string): string => {
  return (
    record.inventory_stock?.material.product_name ??
    record.order?.material.product_name ??
    record.material_usage?.material?.product_name ??
    record.material_usage?.inventory_stock.material.product_name ??
    fallback
  )
}

export const getHistoryRecordUserName = (record: InventoryHistoryListItem, fallback: string): string => {
  return record.performed_by?.label || record.performed_by?.full_name || record.performed_by?.username || fallback
}

export const getHistoryRecordQuantityLabel = (record: InventoryHistoryListItem, fallback: string): string => {
  if (!record.quantity_delta) {
    return fallback
  }

  const unit = record.quantity_unit?.display_name || record.quantity_unit?.unit.label || ''
  return `${formatNumericString(record.quantity_delta)} ${unit}`.trim()
}

export const getHistoryRecordProjectExperimentLabel = (record: InventoryHistoryListItem, fallback: string): string => {
  const project = record.project ?? record.material_usage?.project ?? record.order?.project
  const experiment = record.experiment ?? record.material_usage?.experiment
  const projectName = project?.label || project?.name
  const experimentName = experiment?.label || experiment?.name

  if (projectName && experimentName) {
    return `${projectName} · ${experimentName}`
  }

  return projectName || experimentName || fallback
}

/**
 * Maps an activity to the existing workspace that can show its source record.
 *
 * Returned examples:
 * - `'/inventory/all?preset=all&stock=14'`
 * - `'/inventory/orders?order=8'`
 * - `'/inventory/usages'`
 */
export const getHistoryRecordTargetPath = (record: InventoryHistoryListItem): string | null => {
  if (record.material_usage) {
    return '/inventory/usages'
  }

  if (record.inventory_stock) {
    return `/inventory/all?preset=all&stock=${record.inventory_stock.id}`
  }

  if (record.order) {
    return `/inventory/orders?order=${record.order.id}`
  }

  return null
}
