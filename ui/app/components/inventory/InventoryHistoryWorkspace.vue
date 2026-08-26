<script setup lang="ts">
import { computed, ref } from 'vue'
import InventoryHistoryTable from '~/components/inventory/InventoryHistoryTable.vue'
import { getHistoryRecordTargetPath } from '~/components/inventory/inventory-history.values'
import { useInventoryHistoryPageQuery } from '~/composables/inventory/useInventoryHistoryPageQuery'
import type { InventoryHistoryListItem } from '~/types/inventory'
import { getErrorMessage } from '~/utils/errors'

const pageSize = 20
const { t } = useI18n()
const currentPage = ref(1)
const historyPageQueryParams = computed(() => ({
  page: currentPage.value,
  pageSize,
}))
const historyQuery = useInventoryHistoryPageQuery(historyPageQueryParams)

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
      <p class="inventory-section-title">{{ t('inventory.history_workspace.title') }}</p>
      <p class="text-sm text-slate-600">
        {{ t('inventory.history_workspace.description') }}
      </p>
    </div>

    <UCard :ui="{ root: 'core-card divide-y divide-slate-200/70' }">
      <p v-if="historyQuery.isPending.value" class="text-sm text-slate-600">
        {{ t('inventory.history_workspace.loading') }}
      </p>
      <p v-else-if="errorMessage" class="text-sm text-red-600">
        {{ errorMessage }}
      </p>
      <template v-else>
        <InventoryHistoryTable :records="records" @select-record="openRecord" />

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
