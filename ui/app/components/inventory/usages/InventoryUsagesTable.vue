<script setup lang="ts">
import type { InventoryUsageListItem } from '~/types/inventory'
import { formatDateTime } from '~/utils/dateTime'

type Props = {
  usages: InventoryUsageListItem[]
}

const props = defineProps<Props>()

const toDisplayValue = (value: unknown): string => {
  if (value == null) {
    return '—'
  }

  const text = String(value).trim()
  return text === '' ? '—' : text
}

const usageItemLabel = (usage: InventoryUsageListItem): string => {
  return toDisplayValue(
    usage.material?.label ||
      usage.material?.product_name ||
      usage.inventory_stock?.material?.label ||
      usage.inventory_stock?.material?.product_name,
  )
}
</script>

<template>
  <div class="overflow-x-auto rounded-lg border border-[var(--app-border)] bg-white shadow-sm">
    <table class="min-w-full text-sm">
      <thead class="border-b border-[var(--app-border)] bg-slate-50">
        <tr>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">Product</th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">Used at</th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">Quantity</th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">Unit</th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">Project</th>
          <th class="px-3 py-2 text-left font-semibold text-slate-800">Experiment</th>
        </tr>
      </thead>

      <tbody>
        <tr v-for="usage in props.usages" :key="usage.id" class="border-b border-slate-100 hover:bg-slate-50">
          <td class="px-3 py-2 text-slate-700">
            {{ usageItemLabel(usage) }}
          </td>
          <td class="px-3 py-2 text-slate-700">
            {{ formatDateTime(usage.used_at, { dateStyle: 'medium', timeStyle: 'short' }, '—') }}
          </td>
          <td class="px-3 py-2 text-slate-700">
            {{ toDisplayValue(usage.quantity_used) }}
          </td>
          <td class="px-3 py-2 text-slate-700">
            {{ toDisplayValue(usage.usage_unit.display_name || usage.usage_unit.unit?.label) }}
          </td>
          <td class="px-3 py-2 text-slate-700">
            {{ toDisplayValue(usage.project?.label || usage.project?.name) }}
          </td>
          <td class="px-3 py-2 text-slate-700">
            {{ toDisplayValue(usage.experiment?.label || usage.experiment?.name) }}
          </td>
        </tr>

        <tr v-if="props.usages.length === 0">
          <td colspan="6" class="px-3 py-10 text-center text-sm text-slate-500">No usages found.</td>
        </tr>
      </tbody>
    </table>
  </div>
</template>
