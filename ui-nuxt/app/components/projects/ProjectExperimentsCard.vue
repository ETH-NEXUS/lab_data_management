<script setup lang="ts">
const { t } = useI18n()
</script>

<template>
  <UCard
    class="mx-auto w-[80%]"
    :ui="{
      root: 'border border-white/40 bg-white/30 backdrop-blur-md shadow-sm divide-y divide-white/20',
    }"
  >
    <template #header>
      <p class="font-semibold">{{ t('projects.experiments_card.title') }}</p>
    </template>
    <p class="text-sm text-slate-600">{{ t('projects.experiments_card.placeholder') }}</p>
  </UCard>
</template>
