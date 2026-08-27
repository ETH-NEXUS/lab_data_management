<script setup lang="ts">
import { computed, ref } from 'vue'
import InventoryAwaitingCheckInCard from './InventoryAwaitingCheckInCard.vue'
import InventoryCheckInOutCard from './InventoryCheckInOutCard.vue'
import InventoryDashboardTilePicker from './InventoryDashboardTilePicker.vue'
import InventoryDevicePreviewCard from './InventoryDevicePreviewCard.vue'
import InventoryRecentActivitiesCard from './InventoryRecentActivitiesCard.vue'
import InventoryStockPreviewCard from './InventoryStockPreviewCard.vue'
import InventoryUsagePreviewCard from './InventoryUsagePreviewCard.vue'
import {
  useInventoryStockPageQuery,
  type InventoryStockPageQueryParams,
} from '~/composables/inventory/useInventoryStockPageQuery'
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
const isProjectUsageTileVisible = computed<boolean>(() => visibleDashboardTileKeys.value.has('harvest_project_usages'))
const isExperimentUsageTileVisible = computed<boolean>(() =>
  visibleDashboardTileKeys.value.has('ldm_experiment_usages'),
)
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

const saveDashboardTileKeys = async (tileKeys: string[]): Promise<void> => {
  try {
    await saveTileKeys(tileKeys)
    dashboardTileSaveVersion.value += 1
  } catch {
    // useAPI already shows a readable request error toast.
  }
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

        <InventoryDevicePreviewCard v-if="tile.key === 'device_items'" />

        <InventoryUsagePreviewCard
          v-if="tile.key === 'ldm_experiment_usages'"
          :title="t('inventory.page.actions.recently_linked_ldm_experiments.title')"
          :description="t('inventory.page.actions.recently_linked_ldm_experiments.description')"
          icon="i-heroicons-beaker"
          link-type="experiment"
          :items="recentExperimentUsages"
          :is-loading="recentExperimentUsagesQuery.isPending.value"
          :has-error="Boolean(recentExperimentUsagesQuery.error.value)"
          :error-message="t('inventory.page.actions.recently_linked_ldm_experiments.error')"
        />

        <InventoryUsagePreviewCard
          v-if="tile.key === 'harvest_project_usages'"
          :title="t('inventory.page.actions.recently_linked_harvest_projects.title')"
          :description="t('inventory.page.actions.recently_linked_harvest_projects.description')"
          icon="i-heroicons-link"
          link-type="project"
          :items="recentProjectUsages"
          :is-loading="recentProjectUsagesQuery.isPending.value"
          :has-error="Boolean(recentProjectUsagesQuery.error.value)"
          :error-message="t('inventory.page.actions.recently_linked_harvest_projects.error')"
          footer-text="temporal dev info: the connection to a harvest project happens via material usage - record a material usage and you will see the items here"
        />

        <InventoryAwaitingCheckInCard v-if="tile.key === 'awaiting_check_in'" />
      </template>
    </div>
  </section>
</template>
