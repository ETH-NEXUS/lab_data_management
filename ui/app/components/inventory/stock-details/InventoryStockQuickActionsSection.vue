<script setup lang="ts">
import type { Experiment } from '~/types/experiments'
import type { Project } from '~/types/projects'
import type { InventoryStockListItem } from '~/types/inventory'
import InventoryStockMovePanel from '~/components/inventory/stock-details/InventoryStockMovePanel.vue'
import InventoryStockQuickAdjustPanel from '~/components/inventory/stock-details/InventoryStockQuickAdjustPanel.vue'
import InventoryStockQuickActionsBar from '~/components/inventory/stock-details/InventoryStockQuickActionsBar.vue'
import InventoryStockRecordUsagePanel from '~/components/inventory/stock-details/InventoryStockRecordUsagePanel.vue'
import { useInventoryStockMoveItem } from '~/components/inventory/stock-details/useInventoryStockMoveItem'
import { useInventoryStockQuickAdjust } from '~/components/inventory/stock-details/useInventoryStockQuickAdjust'

type UsageUnitOption = {
  id: number
  label: string
}

type Props = {
  open: boolean
  stock: InventoryStockListItem | null
  isRecordingUsage: boolean
  isSavingUsage: boolean
  selectedProjectId: string
  selectedExperimentId: string
  quantityUsed: string
  selectedUsageUnitId: string
  usageNotes: string
  projects: Project[]
  filteredExperiments: Experiment[]
  usageUnitOptions: UsageUnitOption[]
  stockItemLabel: string
}

const props = defineProps<Props>()

const emit = defineEmits<{
  (e: 'open-usage-mode' | 'cancel-usage-mode' | 'save-usage'): void
  (
    e:
      | 'update:selected-project-id'
      | 'update:selected-experiment-id'
      | 'update:quantity-used'
      | 'update:selected-usage-unit-id'
      | 'update:usage-notes',
    value: string,
  ): void
}>()

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
</script>

<template>
  <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
    <InventoryStockQuickActionsBar
      :is-editing-stock="isEditingStock"
      :is-moving-stock="isMovingStock"
      :is-recording-usage="props.isRecordingUsage"
      @adjust-stock="openEditMode"
      @move-item="openMoveMode"
      @record-usage="emit('open-usage-mode')"
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
      :open="props.isRecordingUsage"
      :is-saving-usage="props.isSavingUsage"
      :stock-item-label="props.stockItemLabel"
      :selected-project-id="props.selectedProjectId"
      :selected-experiment-id="props.selectedExperimentId"
      :quantity-used="props.quantityUsed"
      :selected-usage-unit-id="props.selectedUsageUnitId"
      :usage-notes="props.usageNotes"
      :projects="props.projects"
      :filtered-experiments="props.filteredExperiments"
      :usage-unit-options="props.usageUnitOptions"
      @update:selected-project-id="emit('update:selected-project-id', $event)"
      @update:selected-experiment-id="emit('update:selected-experiment-id', $event)"
      @update:quantity-used="emit('update:quantity-used', $event)"
      @update:selected-usage-unit-id="emit('update:selected-usage-unit-id', $event)"
      @update:usage-notes="emit('update:usage-notes', $event)"
      @cancel="emit('cancel-usage-mode')"
      @save="emit('save-usage')"
    />
  </section>
</template>
