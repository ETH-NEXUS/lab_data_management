<script setup lang="ts">
const { t } = useI18n()
</script>

<template>
  <section class="space-y-3">
    <div class="space-y-1">
      <p class="inventory-section-title">{{ t('inventory.stock_drawer.sections.usage') }}</p>
      <p class="text-sm text-slate-600">Review recorded material usages in a dedicated table workspace.</p>
    </div>
  </section>
</template>
