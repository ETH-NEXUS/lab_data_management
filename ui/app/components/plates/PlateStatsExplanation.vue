<script setup lang="ts">
/**
 * What z' and SSMD below the heatmap mean, with their formulas as utils/plateStats.ts
 * calculates them: from the median and the MAD of the two chosen control well
 * types, e.g. `positive = 'P'`, `negative = 'N'`.
 */
const props = defineProps<{
  positive: string
  negative: string
}>()

const { t } = useI18n()
</script>

<template>
  <div class="mt-2 space-y-1 text-xs text-slate-600">
    <p>{{ t('plates.dynamic.stats.explanation', { positive: props.positive, negative: props.negative }) }}</p>
    <!-- Monospace with spaces, so "-" reads as a minus and not as a dash -->
    <div class="w-fit rounded-md bg-slate-100 px-3 py-2 font-mono text-slate-800 italic">
      <p>{{ t('plates.dynamic.stats.z_prime_formula', { positive: props.positive, negative: props.negative }) }}</p>
      <p>{{ t('plates.dynamic.stats.ssmd_formula', { positive: props.positive, negative: props.negative }) }}</p>
    </div>
    <p>
      <span class="font-medium">{{ t('plates.calculations.reading_caption') }}</span>
      {{ t('plates.dynamic.stats.reading') }}
    </p>
  </div>
</template>
