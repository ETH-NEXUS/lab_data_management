<script setup lang="ts">
import { computed } from 'vue'
import {
  getHistoryRecordActionLabel,
  getHistoryRecordItemName,
  getHistoryRecordRowClass,
  getHistoryRecordTargetPath,
  getHistoryRecordUserName,
} from '~/components/inventory/inventory-history.values'
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

const openRecord = (record: InventoryHistoryListItem): void => {
  const targetPath = getHistoryRecordTargetPath(record)
  if (targetPath) {
    navigateTo(targetPath)
  }
}

const openHistory = (): void => {
  navigateTo('/inventory/history')
}
</script>

<template>
  <UCard :ui="{ root: 'core-card divide-y divide-slate-200/70' }">
    <template #header>
      <div class="flex items-start justify-between gap-3">
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
        <UButton
          variant="ghost"
          color="neutral"
          icon="i-heroicons-arrow-right"
          :aria-label="t('inventory.page.actions.recent_activities.view_all')"
          @click="openHistory"
        />
      </div>
    </template>

    <div class="space-y-2">
      <p v-if="historyQuery.isPending.value" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.loading') }}
      </p>
      <p v-else-if="historyQuery.error.value" class="text-sm text-red-600">
        {{ t('inventory.history_workspace.error') }}
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
          :class="getHistoryRecordRowClass(record)"
          @click="openRecord(record)"
        >
          <span class="absolute -left-[1.05rem] size-3 rounded-full border-2 border-white bg-slate-500" />
          <div class="min-w-0">
            <p class="truncate text-sm font-medium text-slate-800">
              {{ getHistoryRecordActionLabel(record) }} ·
              {{ getHistoryRecordItemName(record, t('inventory.stock_table.values.none')) }}
            </p>
            <p class="truncate text-xs text-slate-600">
              {{ getHistoryRecordUserName(record, 'System') }} ·
              {{ formatDateTime(record.performed_at, { dateStyle: 'medium', timeStyle: 'short' }) }}
            </p>
          </div>
          <UIcon name="i-heroicons-arrow-top-right-on-square" class="size-4 shrink-0 text-slate-500" />
        </button>
      </div>
    </div>
  </UCard>
</template>
