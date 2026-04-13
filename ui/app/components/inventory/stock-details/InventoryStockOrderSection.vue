<script setup lang="ts">
import { computed } from 'vue'
import type { InventoryOrderListItem } from '~/types/inventory'
import { formatDateTime } from '~/utils/dateTime'

type Props = {
  orderEntries: InventoryOrderListItem[]
  isOrdersLoading: boolean
}

const props = defineProps<Props>()

const { t } = useI18n()

const responsibleUsers = computed<string[]>(() => {
  const labels = new Set<string>()

  for (const order of props.orderEntries) {
    const userLabel = order.who_ordered?.label || order.who_ordered?.full_name || order.who_ordered?.username || ''
    if (userLabel !== '') {
      labels.add(userLabel)
    }
  }

  return Array.from(labels)
})
</script>

<template>
  <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
    <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
      {{ t('inventory.stock_drawer.sections.order') }}
    </p>

    <div class="space-y-2 rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
      <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">
        {{ t('inventory.stock_drawer.fields.order_history') }}
      </p>

      <p v-if="props.isOrdersLoading" class="text-sm text-slate-500">
        {{ t('inventory.stock_workspace.loading') }}
      </p>
      <p v-else-if="props.orderEntries.length === 0" class="text-sm text-slate-600">
        {{ t('inventory.stock_drawer.values.no_order_history') }}
      </p>
      <ul v-else class="space-y-1.5 text-sm text-slate-700">
        <li v-for="order in props.orderEntries.slice(0, 8)" :key="order.id">
          <span class="font-medium">
            {{ formatDateTime(order.order_date, { dateStyle: 'medium' }, t('inventory.stock_drawer.values.none')) }}
          </span>
          · {{ order.amount }} {{ order.order_unit.display_name }}
          <span class="text-slate-500"> ({{ order.status_label || order.status }}) </span>
        </li>
      </ul>
    </div>

    <div class="space-y-2 rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
      <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">
        {{ t('inventory.stock_drawer.fields.responsible_user') }}
      </p>
      <p v-if="responsibleUsers.length === 0" class="text-sm text-slate-600">
        {{ t('inventory.stock_drawer.values.no_responsible_user') }}
      </p>
      <ul v-else class="list-disc space-y-1 pl-5 text-sm text-slate-700">
        <li v-for="userLabel in responsibleUsers" :key="userLabel">
          {{ userLabel }}
        </li>
      </ul>
    </div>
  </section>
</template>
