<script setup lang="ts">
import { computed, ref } from 'vue'
import {
  useInventoryStockPageQuery,
  type InventoryStockPageQueryParams,
} from '~/composables/inventory/useInventoryStockPageQuery'
import { useInventoryLookupsQuery } from '~/composables/inventory/useInventoryLookupQuery'
import { useInventoryAwaitingCheckInOrdersQuery } from '~/composables/inventory/useInventoryOrderQuery'
import { useInventoryDashboardTilePreferences } from '~/composables/inventory/useInventoryDashboardTilePreferenceQuery'
import {
  useInventoryRecentExperimentUsagesQuery,
  useInventoryRecentProjectUsagesQuery,
} from '~/composables/inventory/useInventoryUsageQuery'
import { useInventoryStockTablePreferenceStore } from '~/stores/inventory/InventoryStockTablePreferenceStore'
import type { InventoryStockListItem, InventoryStockPreset, InventoryUsageListItem } from '~/types/inventory'

type InventoryPreviewCard = {
  id: string
  title: string
  description: string
  icon: string
  preset: InventoryStockPreset
  items: InventoryStockListItem[]
  isLoading: boolean
  hasError: boolean
}

type VisibleDashboardTile = {
  key: string
  previewCard?: InventoryPreviewCard
}

const { t } = useI18n()
const {
  preferenceQuery: dashboardTilePreferencesQuery,
  isSaving: isSavingDashboardTiles,
  saveTileKeys,
} = useInventoryDashboardTilePreferences()
const dashboardTileSaveVersion = ref(0)
const stockTablePreferenceStore = useInventoryStockTablePreferenceStore()
const dashboardTilePreferences = computed(() => dashboardTilePreferencesQuery.data.value ?? [])
const visibleDashboardTileKeys = computed<Set<string>>(() => {
  return new Set(dashboardTilePreferences.value.filter((tile) => tile.is_visible).map((tile) => tile.key))
})
const isExpiredTileVisible = computed<boolean>(() => visibleDashboardTileKeys.value.has('expired_items'))
const isFavoriteTileVisible = computed<boolean>(() => visibleDashboardTileKeys.value.has('favorite_items'))
const isLowStockTileVisible = computed<boolean>(() => visibleDashboardTileKeys.value.has('low_stock_items'))
const isArchivedTileVisible = computed<boolean>(() => visibleDashboardTileKeys.value.has('archived_items'))
const isDeviceTileVisible = computed<boolean>(() => visibleDashboardTileKeys.value.has('device_items'))
const isAwaitingCheckInTileVisible = computed<boolean>(() => visibleDashboardTileKeys.value.has('awaiting_check_in'))
const isProjectUsageTileVisible = computed<boolean>(() => visibleDashboardTileKeys.value.has('harvest_project_usages'))
const isExperimentUsageTileVisible = computed<boolean>(() =>
  visibleDashboardTileKeys.value.has('ldm_experiment_usages'),
)
const lookupsQuery = useInventoryLookupsQuery(isDeviceTileVisible)
const awaitingCheckInOrdersQuery = useInventoryAwaitingCheckInOrdersQuery(isAwaitingCheckInTileVisible)
const recentProjectUsagesQuery = useInventoryRecentProjectUsagesQuery(isProjectUsageTileVisible)
const recentExperimentUsagesQuery = useInventoryRecentExperimentUsagesQuery(isExperimentUsageTileVisible)
const expiredStockQueryParams = computed<InventoryStockPageQueryParams>(() => ({
  preset: 'expired',
  page: 1,
  pageSize: 5,
  search: '',
  sorting: [{ id: 'expiryDate', desc: true }],
}))
const expiredStocksQuery = useInventoryStockPageQuery(expiredStockQueryParams, isExpiredTileVisible)
const expiredStocks = computed<InventoryStockListItem[]>(() => expiredStocksQuery.data.value?.results ?? [])
const favoriteStockQueryParams = computed<InventoryStockPageQueryParams>(() => ({
  preset: 'favorite',
  page: 1,
  pageSize: 5,
  search: '',
  sorting: stockTablePreferenceStore.sortingState,
}))
const favoriteStocksQuery = useInventoryStockPageQuery(favoriteStockQueryParams, isFavoriteTileVisible)
const favoriteStocks = computed<InventoryStockListItem[]>(() => favoriteStocksQuery.data.value?.results ?? [])
const lowStockQueryParams = computed<InventoryStockPageQueryParams>(() => ({
  preset: 'low_stock',
  page: 1,
  pageSize: 5,
  search: '',
  sorting: stockTablePreferenceStore.sortingState,
}))
const lowStocksQuery = useInventoryStockPageQuery(lowStockQueryParams, isLowStockTileVisible)
const lowStocks = computed<InventoryStockListItem[]>(() => lowStocksQuery.data.value?.results ?? [])
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
const isDeviceStocksQueryEnabled = computed<boolean>(() => isDeviceTileVisible.value && isDeviceSelected.value)
const deviceStocksQuery = useInventoryStockPageQuery(deviceStockQueryParams, isDeviceStocksQueryEnabled)
const deviceStocks = computed<InventoryStockListItem[]>(() => deviceStocksQuery.data.value?.results ?? [])
const awaitingCheckInOrders = computed(() => awaitingCheckInOrdersQuery.data.value ?? [])
const archivedStockQueryParams = computed<InventoryStockPageQueryParams>(() => ({
  preset: 'archived',
  page: 1,
  pageSize: 5,
  search: '',
  sorting: [{ id: 'archivedAt', desc: true }],
}))
const archivedStocksQuery = useInventoryStockPageQuery(archivedStockQueryParams, isArchivedTileVisible)
const archivedStocks = computed<InventoryStockListItem[]>(() => archivedStocksQuery.data.value?.results ?? [])
const recentProjectUsages = computed<InventoryUsageListItem[]>(() => recentProjectUsagesQuery.data.value ?? [])
const recentExperimentUsages = computed<InventoryUsageListItem[]>(() => recentExperimentUsagesQuery.data.value ?? [])

const openOrder = (orderId: number): void => {
  navigateTo(`/inventory/orders?order=${orderId}`)
}

const openUsages = (): void => {
  navigateTo('/inventory/usages')
}

const saveDashboardTileKeys = async (tileKeys: string[]): Promise<void> => {
  try {
    await saveTileKeys(tileKeys)
    dashboardTileSaveVersion.value += 1
  } catch {
    // useAPI already shows a readable request error toast.
  }
}

const getUsageItemName = (usage: InventoryUsageListItem): string => {
  return usage.material?.product_name ?? usage.inventory_stock.material.product_name
}

const previewCards = computed<InventoryPreviewCard[]>(() => [
  {
    id: 'expired_items',
    title: t('inventory.page.actions.recently_expired_items.title'),
    description: t('inventory.page.actions.expired_items.description'),
    icon: 'i-heroicons-clock',
    preset: 'expired',
    items: expiredStocks.value,
    isLoading: expiredStocksQuery.isPending.value,
    hasError: Boolean(expiredStocksQuery.error.value),
  },
  {
    id: 'favorite_items',
    title: t('inventory.page.actions.favorite_items.title'),
    description: t('inventory.page.actions.favorite_items.description'),
    icon: 'i-heroicons-star',
    preset: 'favorite',
    items: favoriteStocks.value,
    isLoading: favoriteStocksQuery.isPending.value,
    hasError: Boolean(favoriteStocksQuery.error.value),
  },
  {
    id: 'low_stock_items',
    title: t('inventory.page.actions.empty_or_low_stock_items.title'),
    description: t('inventory.page.actions.low_stock_items.description'),
    icon: 'i-heroicons-exclamation-triangle',
    preset: 'low_stock',
    items: lowStocks.value,
    isLoading: lowStocksQuery.isPending.value,
    hasError: Boolean(lowStocksQuery.error.value),
  },
  {
    id: 'archived_items',
    title: t('inventory.page.actions.recently_archived_items.title'),
    description: t('inventory.page.actions.archived_items.description'),
    icon: 'i-heroicons-archive-box',
    preset: 'archived',
    items: archivedStocks.value,
    isLoading: archivedStocksQuery.isPending.value,
    hasError: Boolean(archivedStocksQuery.error.value),
  },
])

const visibleDashboardTiles = computed<VisibleDashboardTile[]>(() => {
  const previewCardsById = new Map(previewCards.value.map((card) => [card.id, card]))

  return dashboardTilePreferences.value
    .filter((tile) => tile.is_visible)
    .sort((leftTile, rightTile) => leftTile.position - rightTile.position)
    .map((tile) => ({
      key: tile.key,
      previewCard: previewCardsById.get(tile.key),
    }))
})
</script>

<template>
  <section class="space-y-3">
    <InventoryDashboardTilePicker
      :tiles="dashboardTilePreferences"
      :is-loading="dashboardTilePreferencesQuery.isPending.value"
      :is-saving="isSavingDashboardTiles"
      :has-error="Boolean(dashboardTilePreferencesQuery.error.value)"
      :save-version="dashboardTileSaveVersion"
      @save="saveDashboardTileKeys"
    />

    <p v-if="dashboardTilePreferencesQuery.isPending.value" class="text-sm text-slate-600">
      {{ t('inventory.dashboard_tiles.loading') }}
    </p>
    <p v-else-if="dashboardTilePreferencesQuery.error.value" class="text-sm text-red-600">
      {{ t('inventory.dashboard_tiles.error') }}
    </p>
    <div v-else class="grid gap-4 sm:grid-cols-2 xl:grid-cols-3">
      <p
        v-if="visibleDashboardTiles.length === 0"
        class="rounded-lg border border-dashed border-slate-200 px-4 py-5 text-sm text-slate-600 sm:col-span-2 xl:col-span-3"
      >
        {{ t('inventory.dashboard_tiles.empty_dashboard') }}
      </p>
      <template v-for="tile in visibleDashboardTiles" :key="tile.key">
        <InventoryStockPreviewCard
          v-if="tile.previewCard"
          :key="tile.key"
          :title="tile.previewCard.title"
          :description="tile.previewCard.description"
          :icon="tile.previewCard.icon"
          :preset="tile.previewCard.preset"
          :items="tile.previewCard.items"
          :is-loading="tile.previewCard.isLoading"
          :has-error="tile.previewCard.hasError"
        />

        <InventoryRecentActivitiesCard v-if="tile.key === 'recent_activities'" />

        <InventoryCheckInOutCard v-if="tile.key === 'check_in_out'" />

        <UCard v-if="tile.key === 'device_items'" :ui="{ root: 'core-card divide-y divide-slate-200/70' }">
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

        <UCard v-if="tile.key === 'ldm_experiment_usages'" :ui="{ root: 'core-card divide-y divide-slate-200/70' }">
          <template #header>
            <div class="flex items-start gap-2">
              <span class="inventory-icon-chip">
                <UIcon name="i-heroicons-beaker" class="size-5" />
              </span>
              <div>
                <p class="text-sm font-semibold text-slate-800">
                  {{ t('inventory.page.actions.recently_linked_ldm_experiments.title') }}
                </p>
                <p class="text-sm text-slate-600">
                  {{ t('inventory.page.actions.recently_linked_ldm_experiments.description') }}
                </p>
              </div>
            </div>
          </template>

          <div class="space-y-2">
            <p v-if="recentExperimentUsagesQuery.isPending.value" class="text-sm text-slate-600">
              {{ t('inventory.stock_workspace.loading') }}
            </p>
            <p v-else-if="recentExperimentUsagesQuery.error.value" class="text-sm text-red-600">
              {{ t('inventory.page.actions.recently_linked_ldm_experiments.error') }}
            </p>
            <p v-else-if="recentExperimentUsages.length === 0" class="text-sm text-slate-600">
              {{ t('inventory.stock_workspace.empty') }}
            </p>
            <template v-else>
              <button
                v-for="usage in recentExperimentUsages"
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
                    {{ usage.experiment?.name ?? t('inventory.stock_table.values.none') }}
                  </p>
                  <UIcon name="i-heroicons-arrow-top-right-on-square" class="h-4 w-4 shrink-0 text-slate-400" />
                </div>
              </button>
            </template>
          </div>
        </UCard>

        <UCard v-if="tile.key === 'harvest_project_usages'" :ui="{ root: 'core-card divide-y divide-slate-200/70' }">
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
            <p v-if="recentProjectUsagesQuery.isPending.value" class="text-sm text-slate-600">
              {{ t('inventory.stock_workspace.loading') }}
            </p>
            <p v-else-if="recentProjectUsagesQuery.error.value" class="text-sm text-red-600">
              {{ t('inventory.page.actions.recently_linked_harvest_projects.error') }}
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
              temporal dev info: the connection to a harvest project happens via material usage - record a material
              usage and you will see the items here
            </p>
          </div>
        </UCard>

        <UCard v-if="tile.key === 'awaiting_check_in'" :ui="{ root: 'core-card divide-y divide-slate-200/70' }">
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
      </template>
    </div>
  </section>
</template>
