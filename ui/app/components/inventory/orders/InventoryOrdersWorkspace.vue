<script setup lang="ts">
import { computed, ref } from 'vue'
import { useQueryClient } from '@tanstack/vue-query'
import InventoryOrderCreateModal from '~/components/inventory/orders/InventoryOrderCreateModal.vue'
import InventoryOrdersTable from '~/components/inventory/orders/InventoryOrdersTable.vue'
import { useInventoryMaterialsQuery } from '~/composables/inventory/useInventoryMaterialQuery'
import { useInventoryMaterialOrderUnitsQuery } from '~/composables/inventory/useInventoryMaterialUnitQuery'
import { useInventoryOrdersQuery } from '~/composables/inventory/useInventoryOrderQuery'
import { useInventoryOrderStore } from '~/stores/inventory/InventoryOrderStore'
import { INVENTORY_ORDERS_QUERY_KEY, type CreateInventoryOrderPayload } from '~/types/inventory'
import { getErrorMessage } from '~/utils/errors'

type OrderUnitOption = {
  id: number
  label: string
}

const { t } = useI18n()
const toast = useToast()
const queryClient = useQueryClient()
const inventoryOrderStore = useInventoryOrderStore()

const ordersQuery = useInventoryOrdersQuery()
const materialsQuery = useInventoryMaterialsQuery()

const isCreateOrderModalOpen = ref(false)
const selectedMaterialId = ref<number>(0)

const selectedMaterialIdRef = computed<number>(() => selectedMaterialId.value)
const selectedMaterialOrderUnitsQuery = useInventoryMaterialOrderUnitsQuery(selectedMaterialIdRef)

/**
 * Builds select options for order units from the dedicated material-units endpoint.
 *
 * Accepted data example:
 * - `[{ id: 8, display_name: 'Box', unit: { label: 'Box' } }]`
 *
 * Returned data example:
 * - `[{ id: 8, label: 'Box' }]`
 */
const orderUnitOptions = computed<OrderUnitOption[]>(() => {
  const materialUnits = selectedMaterialOrderUnitsQuery.data.value ?? []
  return materialUnits.map((unit) => ({
    id: unit.id,
    label: unit.display_name || unit.unit?.label || unit.unit?.name || `Unit #${unit.id}`,
  }))
})

const orderUnitsErrorMessage = computed<string | null>(() => {
  const err = selectedMaterialOrderUnitsQuery.error.value
  if (!err) {
    return null
  }

  return getErrorMessage(err)
})

const ordersErrorMessage = computed<string | null>(() => {
  const err = ordersQuery.error.value
  if (!err) {
    return null
  }
  return getErrorMessage(err)
})

const openCreateOrderModal = (): void => {
  selectedMaterialId.value = 0
  isCreateOrderModalOpen.value = true
}

const closeCreateOrderModal = (): void => {
  isCreateOrderModalOpen.value = false
}

/**
 * Updates currently selected material id for dependent order-unit lookup.
 *
 * Accepted data examples:
 * - `12`
 * - `null`
 */
const onSelectMaterial = (materialId: number | null): void => {
  selectedMaterialId.value = materialId ?? 0
}

/**
 * Creates a new inventory order and refreshes order list query.
 *
 * Accepted payload example:
 * - `{ material_id: 12, order_unit_id: 8, amount: '2', order_date: '2026-04-08T00:00:00Z', status: 'ordered' }`
 */
const registerOrder = async (payload: CreateInventoryOrderPayload): Promise<void> => {
  try {
    await inventoryOrderStore.createOrder(payload)
    await queryClient.invalidateQueries({ queryKey: INVENTORY_ORDERS_QUERY_KEY })
    closeCreateOrderModal()

    toast.add({
      title: 'Order registered',
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: 'Failed to register order',
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  }
}
</script>

<template>
  <section class="space-y-3">
    <div class="flex flex-wrap items-center justify-between gap-3">
      <div class="space-y-1">
        <p class="inventory-section-title">{{ t('inventory.page.actions.order.title') }}</p>
        <p class="text-sm text-slate-600">Create new orders and review existing entries.</p>
      </div>

      <UButton
        color="primary"
        icon="i-heroicons-plus"
        label="Register order"
        :disabled="materialsQuery.isPending.value"
        @click="openCreateOrderModal"
      />
    </div>

    <UCard
      :ui="{
        root: 'core-card divide-y divide-slate-200/70',
      }"
    >
      <p v-if="ordersQuery.isPending.value" class="text-sm text-slate-600">Loading inventory orders...</p>
      <p v-else-if="ordersErrorMessage" class="text-sm text-red-600">
        {{ ordersErrorMessage }}
      </p>
      <InventoryOrdersTable v-else :orders="ordersQuery.data.value ?? []" />
    </UCard>

    <InventoryOrderCreateModal
      :open="isCreateOrderModalOpen"
      :materials="materialsQuery.data.value ?? []"
      :order-unit-options="orderUnitOptions"
      :is-materials-loading="materialsQuery.isPending.value"
      :is-order-units-loading="selectedMaterialId > 0 && selectedMaterialOrderUnitsQuery.isPending.value"
      :order-units-error-message="orderUnitsErrorMessage"
      :is-submitting="inventoryOrderStore.isCreatingOrder"
      @update:open="isCreateOrderModalOpen = $event"
      @select-material="onSelectMaterial"
      @submit="registerOrder"
    />
  </section>
</template>
