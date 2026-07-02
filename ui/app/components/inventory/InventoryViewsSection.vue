<script setup lang="ts">
import { computed } from 'vue'
import InventoryActionCard from '~/components/inventory/InventoryActionCard.vue'
import { useInventoryStocksQuery } from '~/composables/inventory/useInventoryStockQuery'
import type { InventoryStockListItem, InventoryStockPreset } from '~/types/inventory'

type InventoryActionItem = {
  id: string
  title: string
  description: string
  icon: string
}

type InventoryPreviewCard = {
  id: string
  title: string
  description: string
  icon: string
  preset: InventoryStockPreset
  items: InventoryStockListItem[]
}

const emit = defineEmits<{
  (e: 'select-action', actionId: string): void
}>()

const { t } = useI18n()
const stocksQuery = useInventoryStocksQuery()
const stocks = computed<InventoryStockListItem[]>(() => stocksQuery.data.value ?? [])

/**
 * Creates lower-grid actions in the existing card visual style.
 *
 * Returned data example:
 * - `[{ id: 'devices', title: 'Devices', description: 'Manage equipment...', icon: 'i-heroicons-computer-desktop' }]`
 */
const inventoryViews = computed<InventoryActionItem[]>(() => [
  {
    id: 'devices',
    title: t('inventory.page.actions.devices.title'),
    description: t('inventory.page.actions.devices.description'),
    icon: 'i-heroicons-computer-desktop',
  },
  {
    id: 'recent_activities',
    title: t('inventory.page.actions.recent_activities.title'),
    description: t('inventory.page.actions.recent_activities.description'),
    icon: 'i-heroicons-bolt',
  },
  {
    id: 'recently_linked_items',
    title: t('inventory.page.actions.recently_linked_items.title'),
    description: t('inventory.page.actions.recently_linked_items.description'),
    icon: 'i-heroicons-link',
  },
])

const onSelectAction = (actionId: string): void => {
  emit('select-action', actionId)
}

const getPreviewItems = (preset: InventoryStockPreset): InventoryStockListItem[] => {
  if (preset === 'favorite') {
    return stocks.value.filter((stock) => stock.is_favorite).slice(0, 5)
  }

  if (preset === 'low_stock') {
    return stocks.value.filter((stock) => stock.is_low_stock || stock.inventory_status === 'low').slice(0, 5)
  }

  if (preset === 'expired') {
    return stocks.value.filter((stock) => stock.is_expired).slice(0, 5)
  }

  return []
}

const previewCards = computed<InventoryPreviewCard[]>(() => [
  {
    id: 'expired_items',
    title: t('inventory.page.actions.expired_items.title'),
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
    title: t('inventory.page.actions.low_stock_items.title'),
    description: t('inventory.page.actions.low_stock_items.description'),
    icon: 'i-heroicons-exclamation-triangle',
    preset: 'low_stock',
    items: getPreviewItems('low_stock'),
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

  if (preset === 'low_stock') {
    return stock.location_label ?? t('inventory.stock_table.values.unknown_location')
  }

  return stock.location_label ?? t('inventory.stock_table.values.unknown_location')
}
</script>

<template>
  <section class="space-y-3">
    <div class="grid gap-4 sm:grid-cols-2 xl:grid-cols-3">
      <UCard
        v-for="card in previewCards"
        :key="card.id"
        :ui="{ root: 'core-card divide-y divide-slate-200/70' }"
      >
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

            <UButton
              variant="ghost"
              color="neutral"
              icon="i-heroicons-arrow-right"
              @click="openPreset(card.preset)"
            />
          </div>
        </template>

        <div class="space-y-2">
          <p v-if="stocksQuery.isPending.value" class="text-sm text-slate-600">
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
              class="flex w-full items-start justify-between gap-3 rounded-md px-2 py-2 text-left hover:bg-slate-50"
              @click="openStockPreviewItem(card.preset, stock.id)"
            >
              <div class="min-w-0">
                <p class="truncate text-sm font-medium text-slate-800">
                  {{ stock.material.product_name }}
                </p>
                <p class="truncate text-xs text-slate-500">
                  {{ getPreviewMeta(stock, card.preset) }}
                </p>
              </div>
              <UIcon
                name="i-heroicons-arrow-top-right-on-square"
                class="mt-0.5 h-4 w-4 shrink-0 text-slate-400"
              />
            </button>
          </template>
        </div>
      </UCard>

      <InventoryActionCard
        v-for="action in inventoryViews"
        :key="action.id"
        :item="action"
        size="view"
        @select="onSelectAction"
      />
    </div>
  </section>
</template>
