<script setup lang="ts">
type Props = {
  zPrime: number | null
  ssmd: number | null
}

const props = defineProps<Props>()
const { t } = useI18n()

const formatStat = (value: number | null): string => (value === null ? '-' : value.toFixed(2))
</script>

<template>
  <div class="mt-3 flex items-center gap-3 text-sm">
    <span class="font-semibold text-slate-700">{{ t('plates.dynamic.stats.z_prime') }}</span>
    <span class="rounded-sm bg-green-100 px-2 py-1 font-mono text-slate-900">{{ formatStat(props.zPrime) }}</span>
    <span class="font-semibold text-slate-700">{{ t('plates.dynamic.stats.ssmd') }}</span>
    <span class="rounded-sm bg-green-100 px-2 py-1 font-mono text-slate-900">{{ formatStat(props.ssmd) }}</span>
  </div>
</template>
