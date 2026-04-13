<script setup lang="ts">
import InventoryOrdersWorkspace from '~/components/inventory/orders/InventoryOrdersWorkspace.vue'

const { t } = useI18n()

const goBackToInventory = (): void => {
  navigateTo('/inventory')
}
</script>

<template>
  <section class="inventory-shell px-2 pb-12 sm:px-3 lg:px-4">
    <div class="w-full space-y-5">
      <UButton
        variant="ghost"
        color="neutral"
        icon="i-heroicons-arrow-left"
        :label="t('inventory.page.title')"
        @click="goBackToInventory"
      />

      <InventoryOrdersWorkspace />
    </div>
  </section>
</template>
