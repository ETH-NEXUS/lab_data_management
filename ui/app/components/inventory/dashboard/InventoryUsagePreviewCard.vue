<script setup lang="ts">
import type { InventoryUsageListItem } from '~/types/inventory'

type Props = {
  title: string
  description: string
  icon: string
  linkType: 'project' | 'experiment'
  items: InventoryUsageListItem[]
  isLoading: boolean
  hasError: boolean
  errorMessage: string
  footerText?: string
}

const props = defineProps<Props>()
const { t } = useI18n()

const getItemName = (usage: InventoryUsageListItem): string => {
  return usage.material?.product_name ?? usage.inventory_stock.material.product_name
}

const getLinkLabel = (usage: InventoryUsageListItem): string => {
  if (props.linkType === 'project') {
    return usage.project.label || usage.project.name
  }

  return usage.experiment?.name ?? t('inventory.stock_table.values.none')
}

const openUsages = (): void => {
  navigateTo('/inventory/usages')
}
</script>

<template>
  <UCard :ui="{ root: 'core-card divide-y divide-slate-200/70' }">
    <template #header>
      <div class="flex items-start gap-2">
        <span class="inventory-icon-chip">
          <UIcon :name="props.icon" class="size-5" />
        </span>
        <div>
          <p class="text-sm font-semibold text-slate-800">
            {{ props.title }}
          </p>
          <p class="text-sm text-slate-600">
            {{ props.description }}
          </p>
        </div>
      </div>
    </template>

    <div class="space-y-2">
      <p v-if="props.isLoading" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.loading') }}
      </p>
      <p v-else-if="props.hasError" class="text-sm text-red-600">
        {{ props.errorMessage }}
      </p>
      <p v-else-if="props.items.length === 0" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.empty') }}
      </p>
      <template v-else>
        <button
          v-for="usage in props.items"
          :key="usage.id"
          type="button"
          class="grid w-full grid-cols-[minmax(0,1fr)_auto] items-center gap-3 rounded-md px-2 py-2 text-left hover:bg-slate-50"
          @click="openUsages"
        >
          <p class="min-w-0 truncate text-sm font-medium text-slate-800">
            {{ getItemName(usage) }}
          </p>
          <div class="flex min-w-0 items-center gap-2">
            <p class="max-w-28 truncate text-right text-xs text-slate-500">
              {{ getLinkLabel(usage) }}
            </p>
            <UIcon name="i-heroicons-arrow-top-right-on-square" class="h-4 w-4 shrink-0 text-slate-400" />
          </div>
        </button>
      </template>
      <p v-if="props.footerText" class="border-t border-slate-100 pt-2 text-xs text-slate-500">
        {{ props.footerText }}
      </p>
    </div>
  </UCard>
</template>
