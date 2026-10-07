<script setup lang="ts">
import { computed } from 'vue'
import type { Plate } from '~/types/lab'
import { formatSummaryNumber, summarizeDataset } from '~/utils/plateDatasets'

/**
 * Two compact lines of statistics of one measurement at one time point, e.g.
 * "64 wells · median 2.866 · min 1.519 · max 5.101" and "Median by well type: N 3.47 · P 2.83".
 */
const props = defineProps<{
  plate: Plate
  label: string
  timestampIndex: number
}>()

const { t } = useI18n()

const summary = computed(() => summarizeDataset(props.plate, props.label, props.timestampIndex))
</script>

<template>
  <div class="text-xs text-slate-600">
    <p>
      {{
        t('plates.calculations.summary', {
          wells: summary.wells,
          median: formatSummaryNumber(summary.median),
          min: formatSummaryNumber(summary.min),
          max: formatSummaryNumber(summary.max),
        })
      }}
    </p>
    <p v-if="summary.byWellType.length > 0">
      {{ t('plates.calculations.median_by_well_type') }}
      <span v-for="(entry, index) in summary.byWellType" :key="`median-${entry.wellType}`">
        <span v-if="index > 0"> · </span>
        <span class="font-medium">{{ entry.wellType }}</span> {{ formatSummaryNumber(entry.median) }}
      </span>
    </p>
  </div>
</template>
