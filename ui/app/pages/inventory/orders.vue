<script setup lang="ts">
import { computed } from 'vue'
import InventoryOrdersWorkspace from '~/components/inventory/orders/InventoryOrdersWorkspace.vue'

const { t } = useI18n()
const route = useRoute()
const router = useRouter()

const goBackToInventory = (): void => {
  navigateTo('/inventory')
}

const selectedOrderId = computed<number | null>(() => {
  const orderValue = route.query.order
  const normalizedOrderValue = Array.isArray(orderValue) ? orderValue[0] : orderValue

  if (!normalizedOrderValue) {
    return null
  }

  const parsedOrderId = Number.parseInt(normalizedOrderValue, 10)
  return Number.isNaN(parsedOrderId) ? null : parsedOrderId
})

const updateSelectedOrderId = async (orderId: number | null): Promise<void> => {
  const nextQuery = { ...route.query }

  if (orderId === null) {
    delete nextQuery.order
  } else {
    nextQuery.order = String(orderId)
  }

  await router.replace({
    query: nextQuery,
  })
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

      <InventoryOrdersWorkspace :initial-order-id="selectedOrderId" @update:selected-order-id="updateSelectedOrderId" />
    </div>
  </section>
</template>
