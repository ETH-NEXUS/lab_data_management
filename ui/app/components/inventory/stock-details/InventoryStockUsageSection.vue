<script setup lang="ts">
import type { InventoryUsageListItem } from '~/types/inventory'
import { formatDateTime } from '~/utils/dateTime'

type Props = {
  usageEntries: InventoryUsageListItem[]
  linkedExperimentProjectLabels: string[]
  isUsagesLoading: boolean
}

const props = defineProps<Props>()

const { t } = useI18n()
</script>

<template>
  <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
    <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
      {{ t('inventory.stock_drawer.sections.usage') }}
    </p>

    <div class="space-y-2 rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
      <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">
        {{ t('inventory.stock_drawer.fields.recent_usages') }}
      </p>

      <p v-if="props.isUsagesLoading" class="text-sm text-slate-500">
        {{ t('inventory.stock_workspace.loading') }}
      </p>
      <p v-else-if="props.usageEntries.length === 0" class="text-sm text-slate-600">
        {{ t('inventory.stock_drawer.values.no_usage') }}
      </p>
      <ul v-else class="space-y-1.5 text-sm text-slate-700">
        <li v-for="usage in props.usageEntries.slice(0, 8)" :key="usage.id">
          <span class="font-medium">
            {{
              formatDateTime(
                usage.used_at,
                { dateStyle: 'medium', timeStyle: 'short' },
                t('inventory.stock_drawer.values.none'),
              )
            }}
          </span>
          · {{ usage.quantity_used }} {{ usage.usage_unit.display_name }}
          <span class="text-slate-500">
            ({{ usage.project?.label || usage.project?.name || t('inventory.stock_drawer.values.none') }})
          </span>
        </li>
      </ul>
    </div>

    <div class="space-y-2 rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
      <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">
        {{ t('inventory.stock_drawer.fields.linked_experiments_projects') }}
      </p>
      <p v-if="props.linkedExperimentProjectLabels.length === 0" class="text-sm text-slate-600">
        {{ t('inventory.stock_drawer.values.no_linked_items') }}
      </p>
      <ul v-else class="list-disc space-y-1 pl-5 text-sm text-slate-700">
        <li v-for="linkedLabel in props.linkedExperimentProjectLabels" :key="linkedLabel">
          {{ linkedLabel }}
        </li>
      </ul>
    </div>
  </section>
</template>
