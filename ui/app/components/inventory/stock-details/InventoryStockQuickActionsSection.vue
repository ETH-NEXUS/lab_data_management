<script setup lang="ts">
import type { Experiment } from '~/types/experiments'
import type { Project } from '~/types/projects'
import type { InventoryRoom, InventorySector, InventoryStockListItem } from '~/types/inventory'
import InventoryStockMovePanel from '~/components/inventory/stock-details/InventoryStockMovePanel.vue'
import InventoryStockQuickAdjustPanel from '~/components/inventory/stock-details/InventoryStockQuickAdjustPanel.vue'
import InventoryStockQuickActionsBar from '~/components/inventory/stock-details/InventoryStockQuickActionsBar.vue'
import InventoryStockRecordUsagePanel from '~/components/inventory/stock-details/InventoryStockRecordUsagePanel.vue'
import { useInventoryStockQuickAdjust } from '~/components/inventory/stock-details/useInventoryStockQuickAdjust'

type UsageUnitOption = {
  id: number
  label: string
}

type Props = {
  open: boolean
  stock: InventoryStockListItem | null
  isMovingStock: boolean
  isSavingMove: boolean
  selectedRoomId: string
  selectedSectorId: string
  rooms: InventoryRoom[]
  filteredSectors: InventorySector[]
  isLookupsLoading: boolean
  lookupsErrorMessage?: string | null
  isMoveSaveDisabled: boolean
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

const props = withDefaults(defineProps<Props>(), {
  lookupsErrorMessage: null,
})

const emit = defineEmits<{
  (
    e: 'open-move-mode' | 'open-usage-mode' | 'cancel-move-mode' | 'save-move' | 'cancel-usage-mode' | 'save-usage',
  ): void
  (
    e:
      | 'update:selected-room-id'
      | 'update:selected-sector-id'
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
</script>

<template>
  <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
    <InventoryStockQuickActionsBar
      :is-editing-stock="isEditingStock"
      :is-moving-stock="props.isMovingStock"
      :is-recording-usage="props.isRecordingUsage"
      @adjust-stock="openEditMode"
      @move-item="emit('open-move-mode')"
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
      :open="props.isMovingStock"
      :is-saving-move="props.isSavingMove"
      :selected-room-id="props.selectedRoomId"
      :selected-sector-id="props.selectedSectorId"
      :rooms="props.rooms"
      :filtered-sectors="props.filteredSectors"
      :is-lookups-loading="props.isLookupsLoading"
      :lookups-error-message="props.lookupsErrorMessage"
      :is-move-save-disabled="props.isMoveSaveDisabled"
      @update:selected-room-id="emit('update:selected-room-id', $event)"
      @update:selected-sector-id="emit('update:selected-sector-id', $event)"
      @cancel="emit('cancel-move-mode')"
      @save="emit('save-move')"
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
