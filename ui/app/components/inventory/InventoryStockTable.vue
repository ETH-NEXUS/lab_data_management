<script setup lang="ts">
import { computed } from 'vue'
import BaseDataTable from '~/components/tables/BaseDataTable.vue'
import type { TableRow } from '~/components/tables/base-data-table.utils'
import { createInventoryStockTableColumns } from '~/components/inventory/inventory-stock-table.columns'
import type { InventoryStockListItem } from '~/types/inventory'

type Props = {
  stocks: InventoryStockListItem[]
}

const props = defineProps<Props>()

const emit = defineEmits<{
  (e: 'select-stock', stock: InventoryStockListItem): void
}>()

const { t } = useI18n()

const tableColumns = computed(() => createInventoryStockTableColumns(t))

/**
 * Forwards one selected table row as its original stock payload.
 *
 * Input example:
 * - `row = { id: 14, material: { product_name: 'PBS Buffer' }, ... }`
 */
const onRowClick = (row: TableRow): void => {
  emit('select-stock', row as InventoryStockListItem)
}
</script>

<template>
  <BaseDataTable
    :data="props.stocks as unknown as TableRow[]"
    :columns="tableColumns"
    :frozen-column-count="5"
    row-clickable
    :global-filter-placeholder="t('table.general.search_placeholder')"
    @row-click="onRowClick"
  />
</template>
