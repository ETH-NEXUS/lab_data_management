<script setup lang="ts">
import { computed, ref } from 'vue'
import {
  useInventoryStockPageQuery,
  type InventoryStockPageQueryParams,
} from '~/composables/inventory/useInventoryStockPageQuery'
import { useInventoryLookupsQuery } from '~/composables/inventory/useInventoryLookupQuery'
import type { InventoryStockListItem } from '~/types/inventory'

const { t } = useI18n()
const lookupsQuery = useInventoryLookupsQuery()
const selectedDeviceTypeId = ref<string>('')

const deviceTypeOptions = computed(() => {
  return (lookupsQuery.data.value?.deviceTypes ?? []).map((deviceType) => ({
    label: deviceType.label || deviceType.name,
    value: String(deviceType.id),
  }))
})

const selectedDeviceId = computed<number | null>(() => {
  const parsedDeviceId = Number.parseInt(selectedDeviceTypeId.value, 10)
  return Number.isInteger(parsedDeviceId) && parsedDeviceId > 0 ? parsedDeviceId : null
})

const isDeviceSelected = computed<boolean>(() => selectedDeviceId.value !== null)
const deviceStockQueryParams = computed<InventoryStockPageQueryParams>(() => ({
  preset: 'all',
  page: 1,
  pageSize: 5,
  search: '',
  sorting: [],
  deviceTypeId: selectedDeviceId.value,
}))
const deviceStocksQuery = useInventoryStockPageQuery(deviceStockQueryParams, isDeviceSelected)
const deviceStocks = computed<InventoryStockListItem[]>(() => deviceStocksQuery.data.value?.results ?? [])

const openStock = (stockId: number): void => {
  navigateTo(`/inventory/all?preset=all&stock=${stockId}`)
}
</script>

<template>
  <UCard :ui="{ root: 'core-card divide-y divide-slate-200/70' }">
    <template #header>
      <div class="space-y-3">
        <div class="flex items-center gap-2">
          <span class="inventory-icon-chip">
            <UIcon name="i-heroicons-cpu-chip" class="size-5" />
          </span>
          <p class="text-sm font-semibold text-slate-800">
            {{ t('inventory.page.actions.specific_for_device.title') }}
          </p>
        </div>
        <USelect
          v-model="selectedDeviceTypeId"
          :items="deviceTypeOptions"
          value-key="value"
          label-key="label"
          :placeholder="t('inventory.page.actions.specific_for_device.placeholder')"
          class="w-full"
        />
      </div>
    </template>

    <div class="space-y-2">
      <p v-if="lookupsQuery.isPending.value" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.loading') }}
      </p>
      <p v-else-if="lookupsQuery.error.value" class="text-sm text-red-600">
        {{ t('inventory.stock_workspace.error') }}
      </p>
      <p v-else-if="!isDeviceSelected" class="text-sm text-slate-600">
        {{ t('inventory.page.actions.specific_for_device.description') }}
      </p>
      <p v-else-if="deviceStocksQuery.isPending.value" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.loading') }}
      </p>
      <p v-else-if="deviceStocksQuery.error.value" class="text-sm text-red-600">
        {{ t('inventory.stock_workspace.error') }}
      </p>
      <p v-else-if="deviceStocks.length === 0" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.empty') }}
      </p>
      <template v-else>
        <button
          v-for="stock in deviceStocks"
          :key="stock.id"
          type="button"
          class="grid w-full grid-cols-[minmax(0,1fr)_auto] items-center gap-3 rounded-md px-2 py-2 text-left hover:bg-slate-50"
          @click="openStock(stock.id)"
        >
          <p class="min-w-0 truncate text-sm font-medium text-slate-800">
            {{ stock.material.product_name }}
          </p>
          <div class="flex min-w-0 items-center gap-2">
            <p class="max-w-28 truncate text-right text-xs text-slate-500">
              {{ stock.location_label ?? t('inventory.stock_table.values.unknown_location') }}
            </p>
            <UIcon name="i-heroicons-arrow-top-right-on-square" class="h-4 w-4 shrink-0 text-slate-400" />
          </div>
        </button>
      </template>
    </div>
  </UCard>
</template>
