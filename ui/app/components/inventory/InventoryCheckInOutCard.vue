<script setup lang="ts">
import { computed } from 'vue'
import {
  getHistoryRecordItemName,
  getHistoryRecordProjectExperimentLabel,
  getHistoryRecordQuantityLabel,
  getHistoryRecordTargetPath,
  getHistoryRecordUserName,
} from '~/components/inventory/inventory-history.values'
import { useInventoryCheckInOutHistoryQuery } from '~/composables/inventory/useInventoryCheckInOutHistoryQuery'
import type { InventoryHistoryListItem } from '~/types/inventory'
import { formatDateTime } from '~/utils/dateTime'

const { t } = useI18n()
const checkInOutHistoryQuery = useInventoryCheckInOutHistoryQuery()
const checkInOutRecords = computed<InventoryHistoryListItem[]>(() => checkInOutHistoryQuery.data.value ?? [])

const isCheckOut = (record: InventoryHistoryListItem): boolean => record.performed_action === 'usage_created'

const openRecord = (record: InventoryHistoryListItem): void => {
  const targetPath = getHistoryRecordTargetPath(record)
  if (targetPath) {
    navigateTo(targetPath)
  }
}

const openHistory = (): void => {
  navigateTo('/inventory/check-in-out-history')
}
</script>

<template>
  <UCard :ui="{ root: 'core-card divide-y divide-slate-200/70' }">
    <template #header>
      <div class="flex items-start justify-between gap-3">
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
        <UButton
          variant="ghost"
          color="neutral"
          icon="i-heroicons-arrow-right"
          :aria-label="t('inventory.page.actions.recent_check_in_out.view_all')"
          @click="openHistory"
        />
      </div>
    </template>

    <div class="space-y-2">
      <p v-if="checkInOutHistoryQuery.isPending.value" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.loading') }}
      </p>
      <p v-else-if="checkInOutHistoryQuery.error.value" class="text-sm text-red-600">
        {{ t('inventory.history_workspace.error') }}
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
            · {{ getHistoryRecordItemName(record, t('inventory.stock_table.values.none')) }}
          </p>
          <p class="truncate text-xs text-slate-600">
            {{ getHistoryRecordUserName(record, 'System') }} ·
            {{ getHistoryRecordQuantityLabel(record, t('inventory.stock_table.values.none')) }} ·
            {{ getHistoryRecordProjectExperimentLabel(record, t('inventory.stock_table.values.none')) }}
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
