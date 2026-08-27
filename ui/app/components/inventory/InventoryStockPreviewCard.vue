<script setup lang="ts">
import { formatNumericString, getStatusLabel } from '~/components/inventory/inventory-stock-table.values'
import type { InventoryStockListItem, InventoryStockPreset } from '~/types/inventory'
import { formatDateTime } from '~/utils/dateTime'

type Props = {
  title: string
  description: string
  icon: string
  preset: InventoryStockPreset
  items: InventoryStockListItem[]
  isLoading: boolean
  hasError: boolean
}

const props = defineProps<Props>()
const { t } = useI18n()

const openPreset = (): void => {
  navigateTo(`/inventory/all?preset=${props.preset}`)
}

const openStock = (stockId: number): void => {
  navigateTo(`/inventory/all?preset=${props.preset}&stock=${stockId}`)
}

const getPreviewMeta = (stock: InventoryStockListItem): string => {
  if (props.preset === 'expired') {
    return stock.expiry_date ?? t('inventory.stock_table.values.none')
  }

  if (props.preset === 'archived') {
    return formatDateTime(stock.archived_at, { dateStyle: 'medium' }, t('inventory.stock_table.values.none'))
  }

  if (props.preset === 'low_stock') {
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
</script>

<template>
  <UCard :ui="{ root: 'core-card divide-y divide-slate-200/70' }">
    <template #header>
      <div class="flex items-start justify-between gap-3">
        <div class="space-y-1">
          <div class="flex items-center gap-2">
            <span class="inventory-icon-chip">
              <UIcon :name="props.icon" class="size-5" />
            </span>
            <p class="text-sm font-semibold text-slate-800">{{ props.title }}</p>
          </div>
          <p class="text-sm text-slate-600">{{ props.description }}</p>
        </div>

        <UButton variant="ghost" color="neutral" icon="i-heroicons-arrow-right" @click="openPreset" />
      </div>
    </template>

    <div class="space-y-2">
      <p v-if="props.isLoading" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.loading') }}
      </p>
      <p v-else-if="props.hasError" class="text-sm text-red-600">
        {{ t('inventory.stock_workspace.error') }}
      </p>
      <p v-else-if="props.items.length === 0" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.empty') }}
      </p>
      <template v-else>
        <button
          v-for="stock in props.items"
          :key="`${props.preset}-${stock.id}`"
          type="button"
          class="grid w-full grid-cols-[minmax(0,1fr)_auto] items-center gap-3 rounded-md px-2 py-2 text-left hover:bg-slate-50"
          @click="openStock(stock.id)"
        >
          <div class="min-w-0">
            <p class="truncate text-sm font-medium text-slate-800">
              {{ stock.material.product_name }}
            </p>
          </div>
          <div class="flex min-w-0 items-center gap-2">
            <UBadge
              v-if="props.preset === 'low_stock'"
              :color="getInventoryStatusColor(stock)"
              variant="soft"
              size="xs"
            >
              {{ getStatusLabel(t, stock.inventory_status) }}
            </UBadge>
            <p class="max-w-28 truncate text-right text-xs text-slate-500">
              {{ getPreviewMeta(stock) }}
            </p>
            <UIcon name="i-heroicons-arrow-top-right-on-square" class="h-4 w-4 shrink-0 text-slate-400" />
          </div>
        </button>
      </template>
    </div>
  </UCard>
</template>
