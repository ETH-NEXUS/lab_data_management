<script setup lang="ts">
import type { InventoryMaterialDetail, InventoryStockListItem } from '~/types/inventory'
import InventoryStockMovePanel from '~/components/inventory/stock-details/InventoryStockMovePanel.vue'
import InventoryStockQuickAdjustPanel from '~/components/inventory/stock-details/InventoryStockQuickAdjustPanel.vue'
import InventoryStockQuickActionsBar from '~/components/inventory/stock-details/InventoryStockQuickActionsBar.vue'
import InventoryStockRecordUsagePanel from '~/components/inventory/stock-details/InventoryStockRecordUsagePanel.vue'
import { useInventoryStockMoveItem } from '~/components/inventory/stock-details/useInventoryStockMoveItem'
import { useInventoryStockQuickAdjust } from '~/components/inventory/stock-details/useInventoryStockQuickAdjust'
import { useInventoryStockRecordUsage } from '~/components/inventory/stock-details/useInventoryStockRecordUsage'

type Props = {
  open: boolean
  stock: InventoryStockListItem | null
  materialDetail: InventoryMaterialDetail | null
}

const props = defineProps<Props>()

const {
  isEditingStock,
  isSavingStockAdjustment,
  editQuantity,
  editMinimumQuantity,
  editNotes,
  stockUnitLabel,
  openEditMode,
  cancelEditMode,
  saveStockAdjustment,
} = useInventoryStockQuickAdjust(props)

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
  <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
    <InventoryStockQuickActionsBar
      :is-editing-stock="isEditingStock"
      :is-moving-stock="isMovingStock"
      :is-recording-usage="isRecordingUsage"
      @adjust-stock="openEditMode"
      @move-item="openMoveMode"
      @record-usage="openUsageMode"
    />

    <InventoryStockQuickAdjustPanel
      :open="isEditingStock"
      :is-saving-stock-adjustment="isSavingStockAdjustment"
      :edit-quantity="editQuantity"
      :edit-minimum-quantity="editMinimumQuantity"
      :edit-notes="editNotes"
      :stock-unit-label="stockUnitLabel"
      @update:edit-quantity="editQuantity = $event"
      @update:edit-minimum-quantity="editMinimumQuantity = $event"
      @update:edit-notes="editNotes = $event"
      @cancel="cancelEditMode"
      @save="saveStockAdjustment"
    />

    <InventoryStockMovePanel
      :open="isMovingStock"
      :is-saving-move="isSavingMove"
      :selected-room-id="selectedRoomId"
      :selected-sector-id="selectedSectorId"
      :rooms="rooms"
      :filtered-sectors="filteredSectors"
      :is-lookups-loading="isLookupsLoading"
      :lookups-error-message="lookupsErrorMessage"
      :is-move-save-disabled="isMoveSaveDisabled"
      @update:selected-room-id="selectedRoomId = $event"
      @update:selected-sector-id="selectedSectorId = $event"
      @cancel="cancelMoveMode"
      @save="saveMove"
    />

    <InventoryStockRecordUsagePanel
      :open="isRecordingUsage"
      :is-saving-usage="isSavingUsage"
      :stock-item-label="stockItemLabel"
      :selected-project-id="selectedProjectId"
      :selected-experiment-id="selectedExperimentId"
      :quantity-used="quantityUsed"
      :selected-usage-unit-id="selectedUsageUnitId"
      :usage-notes="usageNotes"
      :projects="projects"
      :filtered-experiments="filteredExperiments"
      :usage-unit-options="usageUnitOptions"
      @update:selected-project-id="selectedProjectId = $event"
      @update:selected-experiment-id="selectedExperimentId = $event"
      @update:quantity-used="quantityUsed = $event"
      @update:selected-usage-unit-id="selectedUsageUnitId = $event"
      @update:usage-notes="usageNotes = $event"
      @cancel="cancelUsageMode"
      @save="saveUsage"
    />
  </section>
</template>
