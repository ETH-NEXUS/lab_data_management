<script setup lang="ts">
import { computed } from 'vue'
import type { InventoryOrderListItem } from '~/types/inventory'
import { formatDateTime } from '~/utils/dateTime'

type Props = {
  orders: InventoryOrderListItem[]
  selectedOrderId?: number | null
}

const props = defineProps<Props>()
const emit = defineEmits<{
  (e: 'select-order', orderId: number): void
}>()
const { t } = useI18n()

/**
 * Sorts orders by order date descending for newest-first scanning.
 *
 * Returned data example:
 * - `[{ id: 15, order_date: '2026-04-08T10:20:00Z' }, { id: 12, order_date: '2026-04-01T12:00:00Z' }]`
 */
const sortedOrders = computed<InventoryOrderListItem[]>(() => {
  const items = [...props.orders]
  items.sort((leftOrder, rightOrder) => {
    const leftDate = new Date(leftOrder.order_date).getTime()
    const rightDate = new Date(rightOrder.order_date).getTime()
    return rightDate - leftDate
  })
  return items
})

const toDisplayValue = (value: unknown): string => {
  if (value == null) {
    return '—'
  }

  const normalizedValue = String(value).trim()
  if (normalizedValue === '') {
    return '—'
  }

  return normalizedValue
}

const statusBadgeClass = (status: string): string => {
  if (status === 'product_arrived') {
    return 'border-emerald-300 bg-emerald-100 text-emerald-900'
  }

  if (status === 'tentative') {
    return 'border-amber-300 bg-amber-100 text-amber-900'
  }

  return 'border-blue-300 bg-blue-100 text-blue-900'
}
</script>

<template>
  <div class="overflow-x-auto rounded-lg border border-[var(--app-border)] bg-white shadow-sm">
    <table class="min-w-full text-sm">
      <thead class="border-b border-[var(--app-border)] bg-slate-50">
        <tr>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">
            {{ t('inventory.orders.table.columns.order_number') }}
          </th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">
            {{ t('inventory.orders.table.columns.product') }}
          </th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">
            {{ t('inventory.orders.table.columns.date') }}
          </th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">
            {{ t('inventory.orders.table.columns.amount') }}
          </th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">
            {{ t('inventory.orders.table.columns.unit') }}
          </th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">
            {{ t('inventory.orders.table.columns.status') }}
          </th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">
            {{ t('inventory.orders.table.columns.ordered_by') }}
          </th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">
            {{ t('inventory.orders.table.columns.project') }}
          </th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">
            {{ t('inventory.orders.table.columns.notes') }}
          </th>
        </tr>
      </thead>

      <tbody>
        <tr
          v-for="order in sortedOrders"
          :key="order.id"
          :class="[
            'border-b border-slate-100 transition-colors hover:bg-slate-50',
            props.selectedOrderId === order.id ? 'bg-blue-50/50' : '',
          ]"
        >
          <td class="px-3 py-2 text-slate-700">
            {{ toDisplayValue(order.id) }}
          </td>
          <td class="px-3 py-2 text-slate-700">
            <button
              type="button"
              class="cursor-pointer text-left text-sky-700 hover:text-sky-800 hover:underline"
              @click="emit('select-order', order.id)"
            >
              {{ toDisplayValue(order.material.label || order.material.product_name) }}
            </button>
          </td>
          <td class="px-3 py-2 text-slate-700">
            {{ formatDateTime(order.order_date, { dateStyle: 'medium' }, '—') }}
          </td>

          <td class="px-3 py-2 text-slate-700">
            {{ toDisplayValue(order.amount) }}
          </td>
          <td class="px-3 py-2 text-slate-700">
            {{ toDisplayValue(order.order_unit.display_name || order.order_unit.unit?.label) }}
          </td>
          <td class="px-3 py-2">
            <span
              :class="[
                'inline-flex items-center rounded-full border px-2 py-0.5 text-[11px] font-semibold',
                statusBadgeClass(order.status),
              ]"
            >
              {{ order.status_label || order.status }}
            </span>
          </td>
          <td class="px-3 py-2 text-slate-700">
            {{
              toDisplayValue(order.who_ordered?.label || order.who_ordered?.full_name || order.who_ordered?.username)
            }}
          </td>
          <td class="px-3 py-2 text-slate-700">
            {{ toDisplayValue(order.project?.label || order.project?.name) }}
          </td>
          <td class="px-3 py-2 text-slate-700">
            {{ toDisplayValue(order.notes) }}
          </td>
        </tr>

        <tr v-if="sortedOrders.length === 0">
          <td colspan="9" class="px-3 py-10 text-center text-sm text-slate-500">
            {{ t('inventory.orders.table.empty') }}
          </td>
        </tr>
      </tbody>
    </table>
  </div>
</template>
