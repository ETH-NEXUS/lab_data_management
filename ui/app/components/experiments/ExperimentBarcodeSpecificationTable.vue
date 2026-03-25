<script setup lang="ts">
import { computed } from 'vue'
import type { Barcode, BarcodeSpecification } from '~/types/lab'
import { BARCODE_CSV_COLUMNS, generateBarcodes } from '~/utils/barcode'

const props = defineProps<{
  specification: BarcodeSpecification
}>()

const { t } = useI18n()

/**
 * Label metadata for table columns in a stable order.
 *
 * Returned data example:
 * - `[{ key: 'NorthBarcode', label: 'North barcode' }, { key: 'SouthBarcode', label: 'South barcode' }]`
 */
const barcodeColumns = computed<Array<{ key: keyof Barcode; label: string }>>(() => [
  { key: 'NorthBarcode', label: t('experiments.barcodes.table.columns.north') },
  { key: 'SouthBarcode', label: t('experiments.barcodes.table.columns.south') },
  { key: 'EastBarcode', label: t('experiments.barcodes.table.columns.east') },
  { key: 'WestBarcode', label: t('experiments.barcodes.table.columns.west') },
])

/**
 * Generated barcode table rows for one barcode specification.
 *
 * Returned data example:
 * - `[{ NorthBarcode: 'A001_1$A001_1', SouthBarcode: '', EastBarcode: '', WestBarcode: '' }]`
 */
const rows = computed<Barcode[]>(() =>
  generateBarcodes(props.specification.prefix, props.specification.number_of_plates, props.specification.sides),
)

/**
 * Helper kept for clarity when reading table cells.
 *
 * Returned data examples:
 * - `'A001_1$A001_1'`
 * - `'-'`
 */
const getCellValue = (row: Barcode, columnKey: keyof Barcode): string => {
  const value = row[columnKey]
  return value && value.length > 0 ? value : '-'
}

const columnKeyList = BARCODE_CSV_COLUMNS
</script>

<template>
  <article class="space-y-3 rounded-xl p-4">
    <div class="overflow-x-auto rounded-lg border border-slate-200">
      <table class="min-w-full border-collapse text-xs text-slate-700 sm:text-sm">
        <thead class="bg-slate-100/80">
          <tr>
            <th
              v-for="column in barcodeColumns"
              :key="column.key"
              class="border-b border-slate-200 px-3 py-2 text-left font-semibold whitespace-nowrap"
            >
              {{ column.label }}
            </th>
          </tr>
        </thead>

        <tbody>
          <tr
            v-for="(row, rowIndex) in rows"
            :key="`${props.specification.id}-${rowIndex}`"
            class="odd:bg-white even:bg-slate-50/60"
          >
            <td
              v-for="columnKey in columnKeyList"
              :key="`${props.specification.id}-${rowIndex}-${columnKey}`"
              class="border-b border-slate-100 px-3 py-2 font-mono text-[11px] sm:text-xs"
            >
              {{ getCellValue(row, columnKey) }}
            </td>
          </tr>
        </tbody>
      </table>
    </div>
  </article>
</template>
