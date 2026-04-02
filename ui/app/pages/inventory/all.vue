<script setup lang="ts">
import InventoryStockWorkspace from '~/components/inventory/InventoryStockWorkspace.vue'

const { t } = useI18n()

const goBackToInventory = (): void => {
  navigateTo('/inventory')
}
</script>

<template>
  <section class="inventory-shell px-4 pb-12 sm:px-6 lg:px-8">
    <div class="mx-auto w-full max-w-6xl space-y-5">
      <UButton
        variant="ghost"
        color="neutral"
        icon="i-heroicons-arrow-left"
        :label="t('inventory.page.title')"
        @click="goBackToInventory"
      />

      <InventoryStockWorkspace />
    </div>
  </section>
</template>
