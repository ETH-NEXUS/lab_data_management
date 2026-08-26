<script setup lang="ts">
import { computed, ref } from 'vue'
import {
  formatNumericString,
  getStocksForPreset,
  getStatusLabel,
  sortStocksLikeInventoryTable,
} from '~/components/inventory/inventory-stock-table.values'
import {
  useInventoryStockPageQuery,
  type InventoryStockPageQueryParams,
} from '~/composables/inventory/useInventoryStockPageQuery'
import { useInventoryLookupsQuery } from '~/composables/inventory/useInventoryLookupQuery'
import { useInventoryAwaitingCheckInOrdersQuery } from '~/composables/inventory/useInventoryOrderQuery'
import { useInventoryStocksQuery } from '~/composables/inventory/useInventoryStockQuery'
import { useInventoryUsagesQuery } from '~/composables/inventory/useInventoryUsageQuery'
import { useInventoryStockTablePreferenceStore } from '~/stores/inventory/InventoryStockTablePreferenceStore'
import type { InventoryStockListItem, InventoryStockPreset, InventoryUsageListItem } from '~/types/inventory'
import { formatDateTime } from '~/utils/dateTime'

type InventoryPreviewCard = {
  id: string
  title: string
  description: string
  icon: string
  preset: InventoryStockPreset
  items: InventoryStockListItem[]
}

const { t } = useI18n()
const stocksQuery = useInventoryStocksQuery()
const lookupsQuery = useInventoryLookupsQuery()
const awaitingCheckInOrdersQuery = useInventoryAwaitingCheckInOrdersQuery()
const usagesQuery = useInventoryUsagesQuery()
const stockTablePreferenceStore = useInventoryStockTablePreferenceStore()
const stocks = computed<InventoryStockListItem[]>(() => stocksQuery.data.value ?? [])
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
const deviceStocksQuery = useInventoryStockPageQuery(deviceStockQueryParams)
const deviceStocks = computed<InventoryStockListItem[]>(() => deviceStocksQuery.data.value?.results ?? [])
const awaitingCheckInOrders = computed(() => (awaitingCheckInOrdersQuery.data.value ?? []).slice(0, 5))
const archivedStockQueryParams = computed<InventoryStockPageQueryParams>(() => ({
  preset: 'archived',
  page: 1,
  pageSize: 5,
  search: '',
  sorting: [{ id: 'archivedAt', desc: true }],
}))
const archivedStocksQuery = useInventoryStockPageQuery(archivedStockQueryParams)
const archivedStocks = computed<InventoryStockListItem[]>(() => archivedStocksQuery.data.value?.results ?? [])
const recentProjectUsages = computed<InventoryUsageListItem[]>(() => {
  return [...(usagesQuery.data.value ?? [])]
    .sort((leftUsage, rightUsage) => new Date(rightUsage.used_at).getTime() - new Date(leftUsage.used_at).getTime())
    .slice(0, 5)
})

const openOrder = (orderId: number): void => {
  navigateTo(`/inventory/orders?order=${orderId}`)
}

const openUsages = (): void => {
  navigateTo('/inventory/usages')
}

const getUsageItemName = (usage: InventoryUsageListItem): string => {
  return usage.material?.product_name ?? usage.inventory_stock.material.product_name
}

const getPreviewItems = (preset: InventoryStockPreset): InventoryStockListItem[] => {
  const presetStocks = getStocksForPreset(stocks.value, preset)

  if (preset === 'expired') {
    return [...presetStocks]
      .sort((leftStock, rightStock) => (rightStock.expiry_date ?? '').localeCompare(leftStock.expiry_date ?? ''))
      .slice(0, 5)
  }

  return sortStocksLikeInventoryTable(presetStocks, stockTablePreferenceStore.sortingState, t).slice(0, 5)
}

const previewCards = computed<InventoryPreviewCard[]>(() => [
  {
    id: 'expired_items',
    title: t('inventory.page.actions.recently_expired_items.title'),
    description: t('inventory.page.actions.expired_items.description'),
    icon: 'i-heroicons-clock',
    preset: 'expired',
    items: getPreviewItems('expired'),
  },
  {
    id: 'favorite_items',
    title: t('inventory.page.actions.favorite_items.title'),
    description: t('inventory.page.actions.favorite_items.description'),
    icon: 'i-heroicons-star',
    preset: 'favorite',
    items: getPreviewItems('favorite'),
  },
  {
    id: 'low_stock_items',
    title: t('inventory.page.actions.empty_or_low_stock_items.title'),
    description: t('inventory.page.actions.low_stock_items.description'),
    icon: 'i-heroicons-exclamation-triangle',
    preset: 'low_stock',
    items: getPreviewItems('low_stock'),
  },
  {
    id: 'archived_items',
    title: t('inventory.page.actions.recently_archived_items.title'),
    description: t('inventory.page.actions.archived_items.description'),
    icon: 'i-heroicons-archive-box',
    preset: 'archived',
    items: archivedStocks.value,
  },
])

const openPreset = (preset: InventoryStockPreset): void => {
  navigateTo(`/inventory/all?preset=${preset}`)
}

const openStockPreviewItem = (preset: InventoryStockPreset, stockId: number): void => {
  navigateTo(`/inventory/all?preset=${preset}&stock=${stockId}`)
}

const getPreviewMeta = (stock: InventoryStockListItem, preset: InventoryStockPreset): string => {
  if (preset === 'expired') {
    return stock.expiry_date ?? t('inventory.stock_table.values.none')
  }

  if (preset === 'archived') {
    return formatDateTime(stock.archived_at, { dateStyle: 'medium' }, t('inventory.stock_table.values.none'))
  }

  if (preset === 'low_stock') {
    return `${formatNumericString(stock.quantity)} / ${formatNumericString(stock.minimum_quantity)}`
  }

  return stock.location_label ?? t('inventory.stock_table.values.unknown_location')
}

const getInventoryStatusColor = (stock: InventoryStockListItem): 'error' | 'warning' | 'success' => {
  if (stock.inventory_status === 'out_of_stock') {
    return 'error'
  }

  return stock.inventory_status === 'low' ? 'warning' : 'success'
}

const isPreviewLoading = (preset: InventoryStockPreset): boolean => {
  return preset === 'archived' ? archivedStocksQuery.isPending.value : stocksQuery.isPending.value
}
</script>

<template>
  <section class="space-y-3">
    <div class="grid gap-4 sm:grid-cols-2 xl:grid-cols-3">
      <UCard v-for="card in previewCards" :key="card.id" :ui="{ root: 'core-card divide-y divide-slate-200/70' }">
        <template #header>
          <div class="flex items-start justify-between gap-3">
            <div class="space-y-1">
              <div class="flex items-center gap-2">
                <span class="inventory-icon-chip">
                  <UIcon :name="card.icon" class="size-5" />
                </span>
                <p class="text-sm font-semibold text-slate-800">{{ card.title }}</p>
              </div>
              <p class="text-sm text-slate-600">{{ card.description }}</p>
            </div>

            <UButton variant="ghost" color="neutral" icon="i-heroicons-arrow-right" @click="openPreset(card.preset)" />
          </div>
        </template>

        <div class="space-y-2">
          <p v-if="isPreviewLoading(card.preset)" class="text-sm text-slate-600">
            {{ t('inventory.stock_workspace.loading') }}
          </p>
          <p v-else-if="card.items.length === 0" class="text-sm text-slate-600">
            {{ t('inventory.stock_workspace.empty') }}
          </p>
          <template v-else>
            <button
              v-for="stock in card.items"
              :key="`${card.id}-${stock.id}`"
              type="button"
              class="grid w-full grid-cols-[minmax(0,1fr)_auto] items-center gap-3 rounded-md px-2 py-2 text-left hover:bg-slate-50"
              @click="openStockPreviewItem(card.preset, stock.id)"
            >
              <div class="min-w-0">
                <p class="truncate text-sm font-medium text-slate-800">
                  {{ stock.material.product_name }}
                </p>
              </div>
              <div class="flex min-w-0 items-center gap-2">
                <UBadge
                  v-if="card.preset === 'low_stock'"
                  :color="getInventoryStatusColor(stock)"
                  variant="soft"
                  size="xs"
                >
                  {{ getStatusLabel(t, stock.inventory_status) }}
                </UBadge>
                <p class="max-w-28 truncate text-right text-xs text-slate-500">
                  {{ getPreviewMeta(stock, card.preset) }}
                </p>
                <UIcon name="i-heroicons-arrow-top-right-on-square" class="h-4 w-4 shrink-0 text-slate-400" />
              </div>
            </button>
          </template>
        </div>
      </UCard>

      <InventoryRecentActivitiesCard />

      <InventoryCheckInOutCard />

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
          <p v-else-if="!isDeviceSelected" class="text-sm text-slate-600">
            {{ t('inventory.page.actions.specific_for_device.description') }}
          </p>
          <p v-else-if="deviceStocksQuery.isPending.value" class="text-sm text-slate-600">
            {{ t('inventory.stock_workspace.loading') }}
          </p>
          <p v-else-if="deviceStocks.length === 0" class="text-sm text-slate-600">
            {{ t('inventory.stock_workspace.empty') }}
          </p>
          <template v-else>
            <button
              v-for="stock in deviceStocks"
              :key="`device-${stock.id}`"
              type="button"
              class="grid w-full grid-cols-[minmax(0,1fr)_auto] items-center gap-3 rounded-md px-2 py-2 text-left hover:bg-slate-50"
              @click="openStockPreviewItem('all', stock.id)"
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

      <UCard :ui="{ root: 'core-card divide-y divide-slate-200/70' }">
        <template #header>
          <div class="flex items-start gap-2">
            <span class="inventory-icon-chip">
              <UIcon name="i-heroicons-link" class="size-5" />
            </span>
            <div>
              <p class="text-sm font-semibold text-slate-800">
                {{ t('inventory.page.actions.recently_linked_harvest_projects.title') }}
              </p>
              <p class="text-sm text-slate-600">
                {{ t('inventory.page.actions.recently_linked_harvest_projects.description') }}
              </p>
            </div>
          </div>
        </template>

        <div class="space-y-2">
          <p v-if="usagesQuery.isPending.value" class="text-sm text-slate-600">
            {{ t('inventory.stock_workspace.loading') }}
          </p>
          <p v-else-if="recentProjectUsages.length === 0" class="text-sm text-slate-600">
            {{ t('inventory.stock_workspace.empty') }}
          </p>
          <template v-else>
            <button
              v-for="usage in recentProjectUsages"
              :key="usage.id"
              type="button"
              class="grid w-full grid-cols-[minmax(0,1fr)_auto] items-center gap-3 rounded-md px-2 py-2 text-left hover:bg-slate-50"
              @click="openUsages"
            >
              <p class="min-w-0 truncate text-sm font-medium text-slate-800">
                {{ getUsageItemName(usage) }}
              </p>
              <div class="flex min-w-0 items-center gap-2">
                <p class="max-w-28 truncate text-right text-xs text-slate-500">
                  {{ usage.project.label || usage.project.name }}
                </p>
                <UIcon name="i-heroicons-arrow-top-right-on-square" class="h-4 w-4 shrink-0 text-slate-400" />
              </div>
            </button>
          </template>
          <p class="border-t border-slate-100 pt-2 text-xs text-slate-500">
            temporal dev info: the connection to a harvest project happens via material usage - record a material usage
            and you will see the items here
          </p>
        </div>
      </UCard>

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
          <p v-else-if="awaitingCheckInOrders.length === 0" class="text-sm text-slate-600">
            {{ t('inventory.stock_workspace.empty') }}
          </p>
          <div
            v-for="order in awaitingCheckInOrders"
            :key="order.id"
            class="flex items-center gap-2 rounded-md px-2 py-2"
          >
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
    </div>
  </section>
</template>
