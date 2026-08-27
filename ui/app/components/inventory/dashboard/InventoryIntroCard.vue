<script setup lang="ts">
const { t } = useI18n()
</script>

<template>
  <UCard
    :ui="{
      root: 'core-card inventory-intro-card',
    }"
  >
    <div class="space-y-4 sm:space-y-5">
      <p class="text-xs font-semibold tracking-[0.16em] text-slate-500 uppercase">
        {{ t('inventory.page.eyebrow') }}
      </p>

      <h1 class="text-3xl font-semibold tracking-tight text-slate-900 sm:text-4xl">
        {{ t('inventory.page.title') }}
      </h1>

      <p class="max-w-3xl text-sm leading-7 text-slate-600 sm:text-base">
        {{ t('inventory.page.description') }}
      </p>

      <div class="rounded-xl border border-blue-100/80 bg-blue-50/80 px-4 py-3 text-sm text-slate-700">
        {{ t('inventory.page.helper') }}
      </div>
    </div>
  </UCard>
</template>
