<script setup lang="ts">
import { computed, ref } from 'vue'
import InventoryHistoryTable from './InventoryHistoryTable.vue'
import { getHistoryRecordTargetPath } from './inventory-history.values'
import { useInventoryHistoryPageQuery } from '~/composables/inventory/useInventoryHistoryPageQuery'
import type { InventoryHistoryPageQueryParams } from '~/composables/inventory/inventoryHistoryQuery.utils'
import type { InventoryHistoryListItem } from '~/types/inventory'
import { getErrorMessage } from '~/utils/errors'

type Props = {
  variant?: 'activity' | 'check_in_out'
}

const props = withDefaults(defineProps<Props>(), {
  variant: 'activity',
})

const pageSize = 20
const { t } = useI18n()
const currentPage = ref(1)
const isCheckInOut = computed<boolean>(() => props.variant === 'check_in_out')
const historyPageQueryParams = computed<InventoryHistoryPageQueryParams>(() => ({
  page: currentPage.value,
  pageSize,
  activityGroup: isCheckInOut.value ? 'check_in_out' : undefined,
}))
const historyQuery = useInventoryHistoryPageQuery(historyPageQueryParams)

const workspaceTitle = computed<string>(() => {
  return isCheckInOut.value ? t('inventory.check_in_out_workspace.title') : t('inventory.history_workspace.title')
})
const workspaceDescription = computed<string>(() => {
  return isCheckInOut.value
    ? t('inventory.check_in_out_workspace.description')
    : t('inventory.history_workspace.description')
})
const loadingLabel = computed<string>(() => {
  return isCheckInOut.value ? t('inventory.check_in_out_workspace.loading') : t('inventory.history_workspace.loading')
})
const records = computed<InventoryHistoryListItem[]>(() => historyQuery.data.value?.results ?? [])
const totalCount = computed<number>(() => historyQuery.data.value?.count ?? 0)
const totalPages = computed<number>(() => Math.max(1, Math.ceil(totalCount.value / pageSize)))
const hasPreviousPage = computed<boolean>(() => currentPage.value > 1)
const hasNextPage = computed<boolean>(() => historyQuery.data.value?.next != null)
const errorMessage = computed<string | null>(() => {
  return historyQuery.error.value ? getErrorMessage(historyQuery.error.value) : null
})

const openRecord = (record: InventoryHistoryListItem): void => {
  const targetPath = getHistoryRecordTargetPath(record)
  if (targetPath) {
    navigateTo(targetPath)
  }
}

const goToPreviousPage = (): void => {
  if (hasPreviousPage.value) {
    currentPage.value -= 1
  }
}

const goToNextPage = (): void => {
  if (hasNextPage.value) {
    currentPage.value += 1
  }
}
</script>

<template>
  <section class="space-y-3">
    <div class="space-y-1">
      <p class="inventory-section-title">{{ workspaceTitle }}</p>
      <p class="text-sm text-slate-600">
        {{ workspaceDescription }}
      </p>
    </div>

    <UCard :ui="{ root: 'core-card divide-y divide-slate-200/70' }">
      <p v-if="historyQuery.isPending.value" class="text-sm text-slate-600">
        {{ loadingLabel }}
      </p>
      <p v-else-if="errorMessage" class="text-sm text-red-600">
        {{ errorMessage }}
      </p>
      <template v-else>
        <InventoryHistoryTable :records="records" :mode="props.variant" @select-record="openRecord" />

        <div v-if="totalCount > 0" class="mt-4 flex items-center justify-end gap-2">
          <span class="text-xs text-slate-600">
            {{ t('inventory.history_workspace.pagination.page', { current: currentPage, total: totalPages }) }}
          </span>
          <UButton
            size="xs"
            color="neutral"
            variant="ghost"
            icon="i-heroicons-chevron-left"
            :aria-label="t('inventory.history_workspace.pagination.previous')"
            :disabled="!hasPreviousPage || historyQuery.isFetching.value"
            @click="goToPreviousPage"
          />
          <UButton
            size="xs"
            color="neutral"
            variant="ghost"
            icon="i-heroicons-chevron-right"
            :aria-label="t('inventory.history_workspace.pagination.next')"
            :disabled="!hasNextPage || historyQuery.isFetching.value"
            @click="goToNextPage"
          />
        </div>
      </template>
    </UCard>
  </section>
</template>
