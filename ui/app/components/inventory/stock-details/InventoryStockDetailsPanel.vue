<script setup lang="ts">
import type { InventoryStockDetailsPanelProps } from '~/components/inventory/stock-details/useInventoryStockDetailsPanel'
import InventoryStockDetailsHeader from '~/components/inventory/stock-details/InventoryStockDetailsHeader.vue'
import InventoryStockMaterialIdentitySection from '~/components/inventory/stock-details/InventoryStockMaterialIdentitySection.vue'
import InventoryStockMetadataSection from '~/components/inventory/stock-details/InventoryStockMetadataSection.vue'
import InventoryStockOperationalSection from '~/components/inventory/stock-details/InventoryStockOperationalSection.vue'
import InventoryStockOrderSection from '~/components/inventory/stock-details/InventoryStockOrderSection.vue'
import InventoryStockQuickActionsSection from '~/components/inventory/stock-details/InventoryStockQuickActionsSection.vue'
import InventoryStockSupplierSection from '~/components/inventory/stock-details/InventoryStockSupplierSection.vue'
import InventoryStockUnitConversionSection from '~/components/inventory/stock-details/InventoryStockUnitConversionSection.vue'
import InventoryStockUsageSection from '~/components/inventory/stock-details/InventoryStockUsageSection.vue'
import { useInventoryStockDetailsPanel } from '~/components/inventory/stock-details/useInventoryStockDetailsPanel'
import { useInventoryStockMoveItem } from '~/components/inventory/stock-details/useInventoryStockMoveItem'
import { useInventoryStockRecordUsage } from '~/components/inventory/stock-details/useInventoryStockRecordUsage'

type Props = InventoryStockDetailsPanelProps & {
  isMaterialLoading?: boolean
  isUsagesLoading?: boolean
  isOrdersLoading?: boolean
}

const props = withDefaults(defineProps<Props>(), {
  isMaterialLoading: false,
  isUsagesLoading: false,
  isOrdersLoading: false,
})

const emit = defineEmits<{
  (e: 'close'): void
}>()

const {
  inventoryStatusLabel,
  inventoryStatusColor,
  operationalFields,
  unitConversionFields,
  materialIdentityFields,
  supplierFields,
  sortedUsageEntries,
  linkedExperimentProjectLabels,
  sortedOrderEntries,
  responsibleUsers,
  isMetadataExpanded,
  metadataFields,
  toggleMetadata,
} = useInventoryStockDetailsPanel(props)

const {
  isMovingStock,
  isSavingMove,
  selectedRoomId,
  selectedSectorId,
  rooms,
  filteredSectors,
  isLookupsLoading,
  lookupsErrorMessage,
  isMoveSaveDisabled,
  openMoveMode,
  cancelMoveMode,
  saveMove,
} = useInventoryStockMoveItem(props)

const {
  isRecordingUsage,
  isSavingUsage,
  selectedProjectId,
  selectedExperimentId,
  quantityUsed,
  selectedUsageUnitId,
  usageNotes,
  projects,
  filteredExperiments,
  usageUnitOptions,
  stockItemLabel,
  openUsageMode,
  cancelUsageMode,
  saveUsage,
} = useInventoryStockRecordUsage(props)
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
      <InventoryStockQuickActionsSection
        :open="props.open"
        :stock="props.stock"
        :is-moving-stock="isMovingStock"
        :is-saving-move="isSavingMove"
        :selected-room-id="selectedRoomId"
        :selected-sector-id="selectedSectorId"
        :rooms="rooms"
        :filtered-sectors="filteredSectors"
        :is-lookups-loading="isLookupsLoading"
        :lookups-error-message="lookupsErrorMessage"
        :is-move-save-disabled="isMoveSaveDisabled"
        :is-recording-usage="isRecordingUsage"
        :is-saving-usage="isSavingUsage"
        :selected-project-id="selectedProjectId"
        :selected-experiment-id="selectedExperimentId"
        :quantity-used="quantityUsed"
        :selected-usage-unit-id="selectedUsageUnitId"
        :usage-notes="usageNotes"
        :projects="projects"
        :filtered-experiments="filteredExperiments"
        :usage-unit-options="usageUnitOptions"
        :stock-item-label="stockItemLabel"
        @open-move-mode="openMoveMode"
        @open-usage-mode="openUsageMode"
        @update:selected-room-id="selectedRoomId = $event"
        @update:selected-sector-id="selectedSectorId = $event"
        @cancel-move-mode="cancelMoveMode"
        @save-move="saveMove"
        @update:selected-project-id="selectedProjectId = $event"
        @update:selected-experiment-id="selectedExperimentId = $event"
        @update:quantity-used="quantityUsed = $event"
        @update:selected-usage-unit-id="selectedUsageUnitId = $event"
        @update:usage-notes="usageNotes = $event"
        @cancel-usage-mode="cancelUsageMode"
        @save-usage="saveUsage"
      />

      <InventoryStockOperationalSection
        :inventory-status-label="inventoryStatusLabel"
        :inventory-status-color="inventoryStatusColor"
        :fields="operationalFields"
      />

      <InventoryStockUnitConversionSection :fields="unitConversionFields" />

      <InventoryStockMaterialIdentitySection
        :is-material-loading="props.isMaterialLoading"
        :fields="materialIdentityFields"
      />

      <InventoryStockSupplierSection :fields="supplierFields" />

      <InventoryStockUsageSection
        :usage-entries="sortedUsageEntries"
        :linked-experiment-project-labels="linkedExperimentProjectLabels"
        :is-usages-loading="props.isUsagesLoading"
      />

      <InventoryStockOrderSection
        :order-entries="sortedOrderEntries"
        :responsible-users="responsibleUsers"
        :is-orders-loading="props.isOrdersLoading"
      />

      <InventoryStockMetadataSection
        :is-expanded="isMetadataExpanded"
        :fields="metadataFields"
        @toggle="toggleMetadata"
      />
    </div>
  </aside>
</template>
