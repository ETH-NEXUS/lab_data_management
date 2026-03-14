<script setup lang="ts">
const { t } = useI18n()
</script>

<template>
  <UCard class="mx-auto w-[80%]">
    <template #header>
      <p class="font-semibold">{{ t('projects.experiments_card.title') }}</p>
    </template>
    <p class="text-sm text-slate-600">{{ t('projects.experiments_card.placeholder') }}</p>
  </UCard>
</template>
