<script setup lang="ts">
import type { Column } from '@tanstack/vue-table'
import type { TableRow } from '~/components/tables/base-data-table.utils'

type BaseDataTableHeaderFilterProps = {
  column: Column<TableRow, unknown>
  filterSearch: Record<string, string>
  getVisibleFilterOptions: (column: Column<TableRow, unknown>) => string[]
  initializeColumnFilterDraft: (column: Column<TableRow, unknown>) => void
  isDraftValueSelected: (columnId: string, option: string) => boolean
  setDraftValueSelection: (columnId: string, option: string, checked: boolean) => void
  applyColumnFilter: (column: Column<TableRow, unknown>, close?: () => void) => void
  clearColumnFilter: (column: Column<TableRow, unknown>, close?: () => void) => void
}

const props = defineProps<BaseDataTableHeaderFilterProps>()

const { t } = useI18n()

/**
 * Stores one search string for one column filter popover.
 *
 * Input examples:
 * - `''`
 * - `'stock'`
 */
const updateFilterSearch = (value: string | number): void => {
  props.filterSearch[props.column.id] = String(value)
}
</script>

<template>
  <UPopover v-if="props.column.getCanFilter()" :content="{ side: 'bottom', align: 'end' }">
    <UButton
      variant="ghost"
      color="neutral"
      size="xs"
      icon="i-heroicons-funnel"
      :class="props.column.getFilterValue() ? 'text-blue-700' : 'text-slate-500'"
      @click="props.initializeColumnFilterDraft(props.column)"
    />

    <template #content="{ close }">
      <div class="w-64 space-y-3 p-3">
        <div class="flex flex-col gap-1">
          <UButton
            variant="ghost"
            color="neutral"
            icon="i-heroicons-bars-arrow-up"
            :label="t('table.general.sort_asc')"
            block
            @click="
              () => {
                props.column.toggleSorting(false)
                close()
              }
            "
          />
          <UButton
            variant="ghost"
            color="neutral"
            icon="i-heroicons-bars-arrow-down"
            :label="t('table.general.sort_desc')"
            block
            @click="
              () => {
                props.column.toggleSorting(true)
                close()
              }
            "
          />
          <UButton
            variant="ghost"
            color="neutral"
            icon="i-heroicons-arrows-up-down"
            :label="t('table.general.clear_sort')"
            block
            @click="
              () => {
                props.column.clearSorting()
                close()
              }
            "
          />
        </div>

        <UDivider />

        <UInput
          :model-value="props.filterSearch[props.column.id]"
          :placeholder="t('table.general.filter_options_placeholder')"
          icon="i-heroicons-magnifying-glass"
          size="sm"
          @update:model-value="updateFilterSearch"
        />

        <div class="max-h-44 space-y-1 overflow-y-auto pr-1">
          <label
            v-for="option in props.getVisibleFilterOptions(props.column)"
            :key="`${props.column.id}-${option}`"
            class="flex cursor-pointer items-center gap-2 rounded px-1 py-1 hover:bg-slate-50"
          >
            <UCheckbox
              :model-value="props.isDraftValueSelected(props.column.id, option)"
              @update:model-value="
                (checked) => props.setDraftValueSelection(props.column.id, option, Boolean(checked))
              "
            />
            <span class="truncate text-xs text-slate-700">{{ option }}</span>
          </label>

          <p v-if="props.getVisibleFilterOptions(props.column).length === 0" class="px-1 py-1 text-xs text-slate-500">
            {{ t('table.general.no_filter_options') }}
          </p>
        </div>

        <div class="flex items-center gap-2">
          <UButton
            color="primary"
            size="sm"
            :label="t('table.general.apply_filter')"
            @click="props.applyColumnFilter(props.column, close)"
          />
          <UButton
            variant="ghost"
            color="neutral"
            size="sm"
            :label="t('table.general.clear_filter')"
            @click="props.clearColumnFilter(props.column, close)"
          />
        </div>
      </div>
    </template>
  </UPopover>
</template>
