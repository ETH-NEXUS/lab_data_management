<script setup lang="ts">
import type { Table } from '@tanstack/vue-table'
import type { TableRow } from '~/components/tables/base-data-table.utils'
import { getColumnVisibilityLabel } from '~/components/tables/base-data-table.presentation'

type BaseDataTableToolbarProps = {
  table: Table<TableRow>
  globalFilter: string
  globalFilterPlaceholder: string
  enablePagination: boolean
  totalRowsCount: number
}

const props = defineProps<BaseDataTableToolbarProps>()

const emit = defineEmits<{
  'update:globalFilter': [value: string]
  clearAll: []
}>()

const { t } = useI18n()

/**
 * Emits one updated global filter value to the parent table component.
 *
 * Input examples:
 * - `''`
 * - `'tips'`
 */
const updateGlobalFilter = (value: string | number): void => {
  emit('update:globalFilter', String(value))
}

const clearAllFilters = (): void => {
  emit('clearAll')
}
</script>

<template>
  <div class="flex flex-col gap-2 sm:flex-row sm:items-center sm:justify-between">
    <UInput
      :model-value="props.globalFilter"
      :placeholder="props.globalFilterPlaceholder || t('table.general.search_placeholder')"
      icon="i-heroicons-magnifying-glass"
      class="w-full sm:max-w-sm"
      @update:model-value="updateGlobalFilter"
    />

    <div class="flex items-center gap-2">
      <UBadge color="neutral" variant="soft" class="font-medium">
        {{ t('table.general.rows', { count: props.totalRowsCount }) }}
      </UBadge>

      <UPopover :content="{ side: 'bottom', align: 'end' }">
        <UButton
          variant="ghost"
          color="neutral"
          icon="i-heroicons-view-columns"
          :label="t('table.general.columns')"
        />

        <template #content>
          <div class="w-64 space-y-3 p-3">
            <p class="text-xs font-semibold tracking-[0.06em] text-slate-600 uppercase">
              {{ t('table.general.columns') }}
            </p>

            <div class="max-h-56 space-y-1 overflow-y-auto pr-1">
              <label
                v-for="column in props.table.getAllLeafColumns()"
                :key="`visibility-${column.id}`"
                class="flex cursor-pointer items-center gap-2 rounded px-1 py-1 hover:bg-slate-50"
              >
                <UCheckbox
                  :model-value="column.getIsVisible()"
                  @update:model-value="(checked) => column.toggleVisibility(Boolean(checked))"
                />
                <span class="truncate text-xs text-slate-700">{{ getColumnVisibilityLabel(column) }}</span>
              </label>
            </div>
          </div>
        </template>
      </UPopover>

      <div v-if="props.enablePagination && props.totalRowsCount > 0" class="flex items-center gap-1">
        <span class="text-xs text-slate-600">
          {{ props.table.getState().pagination.pageIndex + 1 }} / {{ props.table.getPageCount() }}
        </span>
        <UButton
          size="xs"
          color="neutral"
          variant="ghost"
          icon="i-heroicons-chevron-left"
          :disabled="!props.table.getCanPreviousPage()"
          @click="props.table.previousPage()"
        />
        <UButton
          size="xs"
          color="neutral"
          variant="ghost"
          icon="i-heroicons-chevron-right"
          :disabled="!props.table.getCanNextPage()"
          @click="props.table.nextPage()"
        />
      </div>

      <UButton
        variant="ghost"
        color="neutral"
        icon="i-heroicons-x-mark"
        :label="t('table.general.clear_all')"
        @click="clearAllFilters"
      />
    </div>
  </div>
</template>
