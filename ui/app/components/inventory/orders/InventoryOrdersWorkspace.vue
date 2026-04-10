<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import { useQueryClient } from '@tanstack/vue-query'
import InventoryOrderCreateModal from '~/components/inventory/orders/InventoryOrderCreateModal.vue'
import InventoryOrdersTable from '~/components/inventory/orders/InventoryOrdersTable.vue'
import { useInventoryMaterialsQuery } from '~/composables/inventory/useInventoryMaterialQuery'
import { useInventoryOrdersQuery } from '~/composables/inventory/useInventoryOrderQuery'
import { useProjectsQuery } from '~/composables/useProjectsQuery'
import { useInventoryOrderStore } from '~/stores/inventory/InventoryOrderStore'
import {
  type CreateInventoryOrderPayload,
  INVENTORY_ORDERS_QUERY_KEY,
  type InventoryOrderListItem,
} from '~/types/inventory'
import { formatDateTime } from '~/utils/dateTime'
import { getErrorMessage } from '~/utils/errors'

const { t } = useI18n()
const toast = useToast()
const queryClient = useQueryClient()
const inventoryOrderStore = useInventoryOrderStore()

const ordersQuery = useInventoryOrdersQuery()
const materialsQuery = useInventoryMaterialsQuery()
const projectsQuery = useProjectsQuery()

const isCreateOrderModalOpen = ref(false)
const selectedOrderId = ref<number | null>(null)
const selectedOrderStatus = ref('')
const selectedProjectId = ref('')

const orderStatusOptions = [
  { value: 'ordered', label: 'Ordered' },
  { value: 'tentative', label: 'Tentative' },
  { value: 'product_arrived', label: 'Product arrived' },
]

const orders = computed<InventoryOrderListItem[]>(() => ordersQuery.data.value ?? [])

const selectedOrder = computed<InventoryOrderListItem | null>(() => {
  if (!selectedOrderId.value) return null
  return orders.value.find((order) => order.id === selectedOrderId.value) ?? null
})

const hasSelectedOrder = computed<boolean>(() => selectedOrder.value !== null)

const projectOptions = computed(() => {
  const projects = [...(projectsQuery.data.value ?? [])]
  projects.sort((leftProject, rightProject) => leftProject.name.localeCompare(rightProject.name))
  return projects
})

const projectsErrorMessage = computed<string | null>(() => {
  const err = projectsQuery.error.value
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
  isCreateOrderModalOpen.value = true
}

const closeCreateOrderModal = (): void => {
  isCreateOrderModalOpen.value = false
}

const toDisplayValue = (value: unknown): string => {
  if (value == null) {
    return '—'
  }

  const text = String(value).trim()
  return text === '' ? '—' : text
}

const openOrderDetails = (orderId: number): void => {
  selectedOrderId.value = orderId
}

const closeOrderDetails = (): void => {
  selectedOrderId.value = null
}

const isUpdateStatusDisabled = computed<boolean>(() => {
  if (!selectedOrder.value) return true
  if (inventoryOrderStore.isUpdatingOrder) return true
  if (ordersQuery.isFetching.value) return true
  if (selectedOrderStatus.value.trim() === '') return true
  return selectedOrderStatus.value === selectedOrder.value.status
})

const isUpdateProjectDisabled = computed<boolean>(() => {
  const order = selectedOrder.value
  if (!order) return true
  if (inventoryOrderStore.isUpdatingOrder) return true
  if (projectsQuery.isFetching.value) return true

  const projectIdText = selectedProjectId.value.trim()
  const nextProjectId = projectIdText === '' ? null : Number.parseInt(projectIdText, 10)
  if (nextProjectId !== null && (!Number.isInteger(nextProjectId) || nextProjectId <= 0)) {
    return true
  }

  const currentProjectId = order.project?.id ?? null
  return nextProjectId === currentProjectId
})

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

const updateSelectedOrderStatus = async (): Promise<void> => {
  const order = selectedOrder.value
  if (!order || selectedOrderStatus.value.trim() === '') {
    return
  }

  try {
    await inventoryOrderStore.updateOrder(order.id, { status: selectedOrderStatus.value })
    await queryClient.invalidateQueries({ queryKey: INVENTORY_ORDERS_QUERY_KEY })
    toast.add({
      title: 'Order status updated',
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: 'Failed to update order status',
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  }
}

const updateSelectedOrderProject = async (): Promise<void> => {
  const order = selectedOrder.value
  if (!order) {
    return
  }

  const projectIdText = selectedProjectId.value.trim()
  const nextProjectId = projectIdText === '' ? null : Number.parseInt(projectIdText, 10)

  if (nextProjectId !== null && (!Number.isInteger(nextProjectId) || nextProjectId <= 0)) {
    return
  }

  try {
    await inventoryOrderStore.updateOrder(order.id, { project_id: nextProjectId })
    await queryClient.invalidateQueries({ queryKey: INVENTORY_ORDERS_QUERY_KEY })
    toast.add({
      title: 'Order project updated',
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: 'Failed to update order project',
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  }
}

watch(
  () => selectedOrder.value,
  (order) => {
    selectedOrderStatus.value = order?.status ?? ''
    selectedProjectId.value = order?.project?.id ? String(order.project.id) : ''
  },
)

watch(
  () => selectedOrderId.value,
  (orderId) => {
    if (!orderId) return
    const hasMatch = orders.value.some((order) => order.id === orderId)
    if (!hasMatch) {
      selectedOrderId.value = null
    }
  },
)
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

    <div
      :class="[
        'grid items-start gap-4',
        hasSelectedOrder ? 'xl:grid-cols-[minmax(0,1.45fr)_minmax(0,1fr)]' : 'grid-cols-1',
      ]"
    >
      <UCard
        class="min-w-0"
        :ui="{
          root: 'core-card divide-y divide-slate-200/70',
        }"
      >
        <p v-if="ordersQuery.isPending.value" class="text-sm text-slate-600">Loading inventory orders...</p>
        <p v-else-if="ordersErrorMessage" class="text-sm text-red-600">
          {{ ordersErrorMessage }}
        </p>
        <InventoryOrdersTable
          v-else
          :orders="orders"
          :selected-order-id="selectedOrderId"
          @select-order="openOrderDetails"
        />
      </UCard>

      <aside
        v-if="hasSelectedOrder && selectedOrder"
        class="flex h-full min-h-[28rem] w-full flex-col overflow-hidden rounded-xl border border-[var(--app-border)] bg-[var(--app-surface)] shadow-[0_16px_44px_rgba(15,23,42,0.12)]"
      >
        <header class="flex items-start justify-between gap-3 border-b border-[var(--app-border)] px-5 py-4">
          <div class="min-w-0 space-y-1">
            <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">Order details</p>
            <h2 class="truncate text-lg font-semibold text-slate-900">
              {{ toDisplayValue(selectedOrder.material.label || selectedOrder.material.product_name) }}
            </h2>
            <p class="text-xs text-slate-600">Order #{{ selectedOrder.id }}</p>
          </div>

          <UButton
            size="xs"
            variant="ghost"
            color="neutral"
            icon="i-heroicons-x-mark"
            label="Close"
            @click="closeOrderDetails"
          />
        </header>

        <div class="space-y-4 overflow-y-auto px-5 py-4">
          <section class="grid gap-2 sm:grid-cols-2">
            <div class="rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
              <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">Date</p>
              <p class="mt-1 text-sm text-slate-800">
                {{ formatDateTime(selectedOrder.order_date, { dateStyle: 'medium' }, '—') }}
              </p>
            </div>
            <div class="rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
              <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">Amount</p>
              <p class="mt-1 text-sm text-slate-800">
                {{ selectedOrder.amount }} {{ selectedOrder.order_unit.display_name }}
              </p>
            </div>
            <div class="rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
              <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">Ordered by</p>
              <p class="mt-1 text-sm text-slate-800">
                {{
                  toDisplayValue(
                    selectedOrder.who_ordered?.label ||
                      selectedOrder.who_ordered?.full_name ||
                      selectedOrder.who_ordered?.username,
                  )
                }}
              </p>
            </div>
            <div class="rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
              <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">Project</p>
              <p class="mt-1 text-sm text-slate-800">
                {{ toDisplayValue(selectedOrder.project?.label || selectedOrder.project?.name) }}
              </p>
            </div>
            <div class="rounded-lg border border-slate-200 bg-slate-50 px-3 py-2 sm:col-span-2">
              <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">Notes</p>
              <p class="mt-1 text-sm text-slate-800">{{ toDisplayValue(selectedOrder.notes) }}</p>
            </div>
          </section>

          <section class="space-y-2 rounded-xl border border-slate-200 bg-white p-4">
            <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">Update status</p>
            <select
              v-model="selectedOrderStatus"
              class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
              :disabled="inventoryOrderStore.isUpdatingOrder"
            >
              <option v-for="statusOption in orderStatusOptions" :key="statusOption.value" :value="statusOption.value">
                {{ statusOption.label }}
              </option>
            </select>
            <UButton
              color="primary"
              label="Save status"
              :loading="inventoryOrderStore.isUpdatingOrder"
              :disabled="isUpdateStatusDisabled"
              @click="updateSelectedOrderStatus"
            />
          </section>

          <section class="space-y-2 rounded-xl border border-slate-200 bg-white p-4">
            <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">Bind project</p>
            <select
              v-model="selectedProjectId"
              class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
              :disabled="inventoryOrderStore.isUpdatingOrder || projectsQuery.isPending.value"
            >
              <option value="">{{ projectsQuery.isPending.value ? 'Loading projects...' : 'No project' }}</option>
              <option v-for="project in projectOptions" :key="project.id" :value="String(project.id)">
                {{ project.name }}
              </option>
            </select>
            <p v-if="projectsErrorMessage" class="text-xs text-red-600">{{ projectsErrorMessage }}</p>
            <UButton
              color="primary"
              label="Save project"
              :loading="inventoryOrderStore.isUpdatingOrder"
              :disabled="isUpdateProjectDisabled"
              @click="updateSelectedOrderProject"
            />
          </section>

          <section class="grid gap-2 sm:grid-cols-2">
            <div class="rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
              <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">Created</p>
              <p class="mt-1 text-sm text-slate-800">
                {{ formatDateTime(selectedOrder.created_at, { dateStyle: 'medium', timeStyle: 'short' }, '—') }}
              </p>
            </div>
            <div class="rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
              <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">Updated</p>
              <p class="mt-1 text-sm text-slate-800">
                {{ formatDateTime(selectedOrder.updated_at, { dateStyle: 'medium', timeStyle: 'short' }, '—') }}
              </p>
            </div>
          </section>
        </div>
      </aside>
    </div>

    <InventoryOrderCreateModal
      :open="isCreateOrderModalOpen"
      :materials="materialsQuery.data.value ?? []"
      :is-materials-loading="materialsQuery.isPending.value"
      :is-submitting="inventoryOrderStore.isCreatingOrder"
      @update:open="isCreateOrderModalOpen = $event"
      @submit="registerOrder"
    />
  </section>
</template>
