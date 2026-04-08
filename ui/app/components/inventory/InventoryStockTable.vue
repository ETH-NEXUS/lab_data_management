<script setup lang="ts">
import { computed } from 'vue'
import BaseDataTable from '~/components/tables/BaseDataTable.vue'
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

const onSelectStock = (stock: InventoryStockListItem): void => {
  emit('select-stock', stock)
}

/**
 * Builds columns with a product-name click action that opens stock details.
 *
 * Input callback example:
 * - `(stock) => emit('select-stock', stock)`
 */
const tableColumns = computed(() => createInventoryStockTableColumns(t, onSelectStock))
</script>

<template>
  <BaseDataTable
    :data="props.stocks as unknown as Record<string, unknown>[]"
    :columns="tableColumns"
    :frozen-column-count="5"
    enable-pagination
    :page-size="50"
    :global-filter-placeholder="t('table.general.search_placeholder')"
  />
</template>
