<script setup lang="ts">
import { computed } from 'vue'
import { useInventoryHistoryQuery } from '~/composables/inventory/useInventoryHistoryQuery'
import type { InventoryHistoryListItem } from '~/types/inventory'
import { formatDateTime } from '~/utils/dateTime'

const { t } = useI18n()
const historyQuery = useInventoryHistoryQuery()
const thirtyDaysAgo = Date.now() - 30 * 24 * 60 * 60 * 1000

const recentActivities = computed<InventoryHistoryListItem[]>(() => {
  return (historyQuery.data.value ?? []).filter((record) => {
    return new Date(record.performed_at).getTime() >= thirtyDaysAgo
  })
})

const additionActions = new Set(['stock_created', 'order_created', 'stock_restored'])
const deletionActions = new Set(['stock_deleted', 'order_deleted', 'usage_deleted', 'stock_archived'])

const getActivityRowClass = (record: InventoryHistoryListItem): string => {
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

const getActionLabel = (record: InventoryHistoryListItem): string => {
  return record.performed_action
    .split('_')
    .map((word) => `${word.charAt(0).toUpperCase()}${word.slice(1)}`)
    .join(' ')
}

const getItemName = (record: InventoryHistoryListItem): string => {
  return (
    record.inventory_stock?.material.product_name ??
    record.order?.material.product_name ??
    record.material_usage?.material?.product_name ??
    record.material_usage?.inventory_stock.material.product_name ??
    t('inventory.stock_table.values.none')
  )
}

const getUserName = (record: InventoryHistoryListItem): string => {
  return record.performed_by?.label || record.performed_by?.full_name || record.performed_by?.username || 'System'
}

const openRecord = (record: InventoryHistoryListItem): void => {
  if (record.inventory_stock) {
    navigateTo(`/inventory/all?preset=all&stock=${record.inventory_stock.id}`)
    return
  }

  if (record.order) {
    navigateTo(`/inventory/orders?order=${record.order.id}`)
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
          <UIcon name="i-heroicons-bolt" class="size-5" />
        </span>
        <div>
          <p class="text-sm font-semibold text-slate-800">
            {{ t('inventory.page.actions.recent_activities.title') }}
          </p>
          <p class="text-sm text-slate-600">
            {{ t('inventory.page.actions.recent_activities.description') }}
          </p>
        </div>
      </div>
    </template>

    <div class="space-y-2">
      <p v-if="historyQuery.isPending.value" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.loading') }}
      </p>
      <p v-else-if="recentActivities.length === 0" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.empty') }}
      </p>
      <div v-else class="space-y-2 border-l-2 border-slate-200 pl-3">
        <button
          v-for="record in recentActivities"
          :key="record.id"
          type="button"
          class="relative grid w-full grid-cols-[minmax(0,1fr)_auto] items-center gap-3 rounded-md px-3 py-2 text-left transition-colors"
          :class="getActivityRowClass(record)"
          @click="openRecord(record)"
        >
          <span class="absolute -left-[1.05rem] size-3 rounded-full border-2 border-white bg-slate-500" />
          <div class="min-w-0">
            <p class="truncate text-sm font-medium text-slate-800">
              {{ getActionLabel(record) }} · {{ getItemName(record) }}
            </p>
            <p class="truncate text-xs text-slate-600">
              {{ getUserName(record) }} ·
              {{ formatDateTime(record.performed_at, { dateStyle: 'medium', timeStyle: 'short' }) }}
            </p>
          </div>
          <UIcon name="i-heroicons-arrow-top-right-on-square" class="size-4 shrink-0 text-slate-500" />
        </button>
      </div>
    </div>
  </UCard>
</template>
