<script setup lang="ts">
import { computed } from 'vue'
import { formatNumericString } from '~/components/inventory/inventory-stock-table.values'
import { useInventoryCheckInOutHistoryQuery } from '~/composables/inventory/useInventoryCheckInOutHistoryQuery'
import type { InventoryHistoryListItem } from '~/types/inventory'
import { formatDateTime } from '~/utils/dateTime'

const { t } = useI18n()
const checkInOutHistoryQuery = useInventoryCheckInOutHistoryQuery()
const checkInOutRecords = computed<InventoryHistoryListItem[]>(() => checkInOutHistoryQuery.data.value ?? [])

const isCheckOut = (record: InventoryHistoryListItem): boolean => record.performed_action === 'usage_created'

const getItemName = (record: InventoryHistoryListItem): string => {
  return (
    record.inventory_stock?.material.product_name ??
    record.material_usage?.material?.product_name ??
    record.material_usage?.inventory_stock.material.product_name ??
    t('inventory.stock_table.values.none')
  )
}

const getUserName = (record: InventoryHistoryListItem): string => {
  return record.performed_by?.label || record.performed_by?.full_name || record.performed_by?.username || 'System'
}

const getQuantityLabel = (record: InventoryHistoryListItem): string => {
  const quantity = record.quantity_delta ? formatNumericString(record.quantity_delta) : '—'
  const unit = record.quantity_unit?.display_name || record.quantity_unit?.unit.label || ''
  return `${quantity} ${unit}`.trim()
}

const getProjectExperimentLabel = (record: InventoryHistoryListItem): string => {
  const project = record.project ?? record.material_usage?.project ?? record.order?.project
  const experiment = record.experiment ?? record.material_usage?.experiment
  const projectName = project?.label || project?.name
  const experimentName = experiment?.label || experiment?.name

  if (projectName && experimentName) {
    return `${projectName} · ${experimentName}`
  }

  return projectName || experimentName || '—'
}

const openRecord = (record: InventoryHistoryListItem): void => {
  if (record.inventory_stock) {
    navigateTo(`/inventory/all?preset=all&stock=${record.inventory_stock.id}`)
    return
  }

  if (record.material_usage) {
    navigateTo('/inventory/usages')
  }
}
</script>

<template>
  <UCard :ui="{ root: 'core-card divide-y divide-slate-200/70' }">
    <template #header>
      <div class="flex items-start gap-2">
        <span class="inventory-icon-chip">
          <UIcon name="i-heroicons-arrow-left-right" class="size-5" />
        </span>
        <div>
          <p class="text-sm font-semibold text-slate-800">
            {{ t('inventory.page.actions.recent_check_in_out.title') }}
          </p>
          <p class="text-sm text-slate-600">
            {{ t('inventory.page.actions.recent_check_in_out.description') }}
          </p>
        </div>
      </div>
    </template>

    <div class="space-y-2">
      <p v-if="checkInOutHistoryQuery.isPending.value" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.loading') }}
      </p>
      <p v-else-if="checkInOutRecords.length === 0" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.empty') }}
      </p>
      <button
        v-for="record in checkInOutRecords"
        :key="record.id"
        type="button"
        class="grid w-full grid-cols-[minmax(0,1fr)_auto] items-center gap-3 rounded-md px-3 py-2 text-left transition-colors"
        :class="isCheckOut(record) ? 'bg-[#F6E7C8] hover:bg-[#F1DDB8]' : 'bg-[#DDEBD8] hover:bg-[#D1E3CB]'"
        @click="openRecord(record)"
      >
        <div class="min-w-0">
          <p class="truncate text-sm font-medium text-slate-800">
            {{
              isCheckOut(record)
                ? t('inventory.page.actions.recent_check_in_out.check_out')
                : t('inventory.page.actions.recent_check_in_out.check_in')
            }}
            · {{ getItemName(record) }}
          </p>
          <p class="truncate text-xs text-slate-600">
            {{ getUserName(record) }} · {{ getQuantityLabel(record) }} · {{ getProjectExperimentLabel(record) }}
          </p>
          <p class="text-xs text-slate-500">
            {{ formatDateTime(record.performed_at, { dateStyle: 'medium', timeStyle: 'short' }) }}
          </p>
        </div>
        <UIcon name="i-heroicons-arrow-top-right-on-square" class="size-4 shrink-0 text-slate-500" />
      </button>
    </div>
  </UCard>
</template>
