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
      <p class="text-primary font-semibold">{{ t('messages_page.sections.empty_wells.title') }}</p>
    </template>

    <div class="max-w-[900px] space-y-3 text-sm text-slate-700">
      <p>{{ t('messages_page.sections.empty_wells.paragraph_1') }}</p>
      <p>{{ t('messages_page.sections.empty_wells.paragraph_2') }}</p>
    </div>
  </UCard>
</template>
