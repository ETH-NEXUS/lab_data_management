<script setup lang="ts">
import { computed } from 'vue'
import InventoryUsagesTable from '~/components/inventory/usages/InventoryUsagesTable.vue'
import { useInventoryUsagesQuery } from '~/composables/inventory/useInventoryUsageQuery'
import type { InventoryUsageListItem } from '~/types/inventory'
import { getErrorMessage } from '~/utils/errors'

const { t } = useI18n()
const usagesQuery = useInventoryUsagesQuery()

const usages = computed<InventoryUsageListItem[]>(() => {
  const items = [...(usagesQuery.data.value ?? [])]
  items.sort((leftUsage, rightUsage) => new Date(rightUsage.used_at).getTime() - new Date(leftUsage.used_at).getTime())
  return items
})

const usagesErrorMessage = computed<string | null>(() => {
  const err = usagesQuery.error.value
  if (!err) {
    return null
  }
  return getErrorMessage(err)
})
</script>

<template>
  <section class="space-y-3">
    <div class="space-y-1">
      <p class="inventory-section-title">{{ t('inventory.stock_drawer.sections.usage') }}</p>
      <p class="text-sm text-slate-600">Review recorded material usages in a dedicated table workspace.</p>
    </div>

    <UCard
      :ui="{
        root: 'core-card divide-y divide-slate-200/70',
      }"
    >
      <p v-if="usagesQuery.isPending.value" class="text-sm text-slate-600">
        {{ t('inventory.stock_workspace.loading') }}
      </p>
      <p v-else-if="usagesErrorMessage" class="text-sm text-red-600">
        {{ usagesErrorMessage }}
      </p>
      <InventoryUsagesTable v-else :usages="usages" />
    </UCard>
  </section>
</template>
