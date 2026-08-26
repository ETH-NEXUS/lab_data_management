<script setup lang="ts">
import {
  getHistoryRecordActionLabel,
  getHistoryRecordItemName,
  getHistoryRecordProjectExperimentLabel,
  getHistoryRecordQuantityLabel,
  getHistoryRecordRowClass,
  getHistoryRecordUserName,
} from '~/components/inventory/inventory-history.values'
import type { InventoryHistoryListItem } from '~/types/inventory'
import { formatDateTime } from '~/utils/dateTime'

type Props = {
  records: InventoryHistoryListItem[]
}

const props = defineProps<Props>()
const emit = defineEmits<{
  (e: 'select-record', record: InventoryHistoryListItem): void
}>()
const { t } = useI18n()

const selectRecord = (record: InventoryHistoryListItem): void => {
  emit('select-record', record)
}
</script>

<template>
  <div class="overflow-x-auto rounded-lg border border-[var(--app-border)] bg-white shadow-sm">
    <table class="min-w-full text-sm">
      <thead class="border-b border-[var(--app-border)] bg-slate-50">
        <tr>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">
            {{ t('inventory.history_workspace.table.action') }}
          </th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">
            {{ t('inventory.history_workspace.table.item') }}
          </th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">
            {{ t('inventory.history_workspace.table.user') }}
          </th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">
            {{ t('inventory.history_workspace.table.quantity') }}
          </th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">
            {{ t('inventory.history_workspace.table.context') }}
          </th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">
            {{ t('inventory.history_workspace.table.performed_at') }}
          </th>
          <th>
            <span class="sr-only">{{ t('inventory.history_workspace.table.open') }}</span>
          </th>
        </tr>
      </thead>

      <tbody>
        <tr
          v-for="record in props.records"
          :key="record.id"
          class="cursor-pointer border-b border-white/70 transition-colors"
          :class="getHistoryRecordRowClass(record)"
          @click="selectRecord(record)"
        >
          <td class="px-3 py-2 font-medium whitespace-nowrap text-slate-800">
            {{ getHistoryRecordActionLabel(record) }}
          </td>
          <td class="max-w-52 truncate px-3 py-2 text-slate-700">
            {{ getHistoryRecordItemName(record, t('inventory.stock_table.values.none')) }}
          </td>
          <td class="px-3 py-2 whitespace-nowrap text-slate-700">
            {{ getHistoryRecordUserName(record, 'System') }}
          </td>
          <td class="px-3 py-2 whitespace-nowrap text-slate-700">
            {{ getHistoryRecordQuantityLabel(record, t('inventory.stock_table.values.none')) }}
          </td>
          <td class="max-w-60 truncate px-3 py-2 text-slate-700">
            {{ getHistoryRecordProjectExperimentLabel(record, t('inventory.stock_table.values.none')) }}
          </td>
          <td class="px-3 py-2 whitespace-nowrap text-slate-700">
            {{ formatDateTime(record.performed_at, { dateStyle: 'medium', timeStyle: 'short' }) }}
          </td>
          <td class="px-3 py-2 text-slate-500">
            <UIcon name="i-heroicons-arrow-top-right-on-square" class="size-4" />
          </td>
        </tr>

        <tr v-if="props.records.length === 0">
          <td colspan="7" class="px-3 py-10 text-center text-sm text-slate-500">
            {{ t('inventory.history_workspace.empty') }}
          </td>
        </tr>
      </tbody>
    </table>
  </div>
</template>
