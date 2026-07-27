<script setup lang="ts">
import { computed } from 'vue'
import InventoryStockDetailsHeader from '~/components/inventory/stock-details/InventoryStockDetailsHeader.vue'
import InventoryStockMaterialIdentitySection from '~/components/inventory/stock-details/InventoryStockMaterialIdentitySection.vue'
import InventoryStockMetadataSection from '~/components/inventory/stock-details/InventoryStockMetadataSection.vue'
import InventoryStockOperationalSection from '~/components/inventory/stock-details/InventoryStockOperationalSection.vue'
import InventoryStockOrderSection from '~/components/inventory/stock-details/InventoryStockOrderSection.vue'
import InventoryStockQuickActionsSection from '~/components/inventory/stock-details/InventoryStockQuickActionsSection.vue'
import InventoryStockSupplierSection from '~/components/inventory/stock-details/InventoryStockSupplierSection.vue'
import InventoryStockUnitConversionSection from '~/components/inventory/stock-details/InventoryStockUnitConversionSection.vue'
import InventoryStockUsageSection from '~/components/inventory/stock-details/InventoryStockUsageSection.vue'
import { useInventoryMaterialQuery } from '~/composables/inventory/useInventoryMaterialQuery'
import { useInventoryOrdersQuery } from '~/composables/inventory/useInventoryOrderQuery'
import { useInventoryUsagesQuery } from '~/composables/inventory/useInventoryUsageQuery'
import { useInventoryStockDetailsPanel } from '~/components/inventory/stock-details/useInventoryStockDetailsPanel'
import type { InventoryOrderListItem, InventoryStockListItem, InventoryUsageListItem } from '~/types/inventory'

type Props = {
  open: boolean
  stock: InventoryStockListItem | null
}

const props = defineProps<Props>()

const selectedMaterialId = computed<number>(() => {
  return props.stock?.material?.id ?? 0
})

const emit = defineEmits<{
  (e: 'close'): void
}>()

const selectedMaterialQuery = useInventoryMaterialQuery(selectedMaterialId)
const usagesQuery = useInventoryUsagesQuery()
const ordersQuery = useInventoryOrdersQuery()

const selectedStockUsageEntries = computed<InventoryUsageListItem[]>(() => {
  const stock = props.stock
  if (!stock) {
    return []
  }

  const usages = usagesQuery.data.value ?? []
  return usages.filter((usage) => usage.inventory_stock.id === stock.id)
})

const selectedStockOrderEntries = computed<InventoryOrderListItem[]>(() => {
  const stock = props.stock
  if (!stock) {
    return []
  }

  const orders = ordersQuery.data.value ?? []
  return orders.filter((order) => order.material.id === stock.material.id)
})

const {
  inventoryStatusLabel,
  inventoryStatusColor,
  storageTemperatureLabel,
  operationalFields,
  unitConversionFields,
  materialIdentityFields,
  supplierFields,
  sortedUsageEntries,
  linkedExperimentProjectLabels,
  sortedOrderEntries,
  isMetadataExpanded,
  metadataFields,
  toggleMetadata,
} = useInventoryStockDetailsPanel({
  get open() {
    return props.open
  },
  get stock() {
    return props.stock
  },
  get materialDetail() {
    return selectedMaterialQuery.data.value ?? null
  },
  get usageEntries() {
    return selectedStockUsageEntries.value
  },
  get orderEntries() {
    return selectedStockOrderEntries.value
  },
})
</script>

<template>
  <aside
    v-if="props.open && props.stock"
    class="flex h-full min-h-[36rem] w-full flex-col overflow-hidden rounded-xl border border-[var(--app-border)] bg-[var(--app-surface)] shadow-[0_16px_44px_rgba(15,23,42,0.12)]"
  >
    <InventoryStockDetailsHeader
      :stock-id="props.stock.id"
      :product-name="props.stock.material.product_name"
      @close="emit('close')"
    />

    <div class="space-y-5 overflow-y-auto px-5 py-4">
      <InventoryStockQuickActionsSection :open="props.open" :stock="props.stock" @close="emit('close')" />

      <InventoryStockOperationalSection
        :inventory-status-label="inventoryStatusLabel"
        :inventory-status-color="inventoryStatusColor"
        :storage-temperature-label="storageTemperatureLabel"
        :fields="operationalFields"
      />

      <InventoryStockUnitConversionSection :fields="unitConversionFields" />

      <InventoryStockMaterialIdentitySection
        :is-material-loading="selectedMaterialQuery.isPending.value"
        :fields="materialIdentityFields"
      />

      <InventoryStockSupplierSection :fields="supplierFields" />

      <InventoryStockUsageSection
        :usage-entries="sortedUsageEntries"
        :linked-experiment-project-labels="linkedExperimentProjectLabels"
        :is-usages-loading="usagesQuery.isPending.value"
      />

      <InventoryStockOrderSection
        :source-order="props.stock.source_order"
        :order-entries="sortedOrderEntries"
        :is-orders-loading="ordersQuery.isPending.value"
      />

      <InventoryStockMetadataSection
        :is-expanded="isMetadataExpanded"
        :fields="metadataFields"
        @toggle="toggleMetadata"
      />
    </div>
  </aside>
</template>
