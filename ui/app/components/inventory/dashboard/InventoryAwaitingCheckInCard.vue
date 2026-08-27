<script setup lang="ts">
import { computed } from 'vue'
import { useInventoryAwaitingCheckInOrdersQuery } from '~/composables/inventory/useInventoryOrderQuery'
import type { InventoryOrderListItem } from '~/types/inventory'

const { t } = useI18n()
const awaitingCheckInOrdersQuery = useInventoryAwaitingCheckInOrdersQuery()
const awaitingCheckInOrders = computed<InventoryOrderListItem[]>(() => awaitingCheckInOrdersQuery.data.value ?? [])

const openOrder = (orderId: number): void => {
  navigateTo(`/inventory/orders?order=${orderId}`)
}
</script>

<template>
  <UCard :ui="{ root: 'core-card divide-y divide-slate-200/70' }">
    <template #header>
      <div class="flex items-start gap-2">
        <span class="inventory-icon-chip">
          <UIcon name="i-heroicons-truck" class="size-5" />
        </span>
        <div>
          <p class="text-sm font-semibold text-slate-800">
            {{ t('inventory.page.actions.awaiting_check_in.title') }}
          </p>
          <p class="text-sm text-slate-600">{{ t('inventory.page.actions.awaiting_check_in.description') }}</p>
        </div>
      </div>
    </template>

    <div class="space-y-2">
      <p v-if="awaitingCheckInOrdersQuery.isPending.value" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.loading') }}
      </p>
      <p v-else-if="awaitingCheckInOrdersQuery.error.value" class="text-sm text-red-600">
        {{ t('inventory.page.actions.awaiting_check_in.error') }}
      </p>
      <p v-else-if="awaitingCheckInOrders.length === 0" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.empty') }}
      </p>
      <div v-for="order in awaitingCheckInOrders" :key="order.id" class="flex items-center gap-2 rounded-md px-2 py-2">
        <p class="min-w-0 flex-1 truncate text-sm font-medium text-slate-800">
          {{ order.material.product_name }}
        </p>
        <UButton
          variant="ghost"
          color="neutral"
          icon="i-heroicons-arrow-top-right-on-square"
          :aria-label="t('inventory.page.actions.awaiting_check_in.open_order')"
          :title="t('inventory.page.actions.awaiting_check_in.open_order')"
          @click="openOrder(order.id)"
        />
      </div>
    </div>
  </UCard>
</template>
