<script setup lang="ts">
import { computed } from 'vue'
import { useQueryClient } from '@tanstack/vue-query'
import type { InventoryStockListItem } from '~/types/inventory'
import InventoryStockMovePanel from '~/components/inventory/stock-details/InventoryStockMovePanel.vue'
import InventoryStockQuickAdjustPanel from '~/components/inventory/stock-details/InventoryStockQuickAdjustPanel.vue'
import InventoryStockQuickActionsBar from '~/components/inventory/stock-details/InventoryStockQuickActionsBar.vue'
import InventoryStockRecordUsagePanel from '~/components/inventory/stock-details/InventoryStockRecordUsagePanel.vue'
import { useInventoryStockStore } from '~/stores/inventory/InventoryStock'
import { useInventoryMaterialQuery } from '~/composables/inventory/useInventoryMaterialQuery'
import { useInventoryStockMoveItem } from '~/components/inventory/stock-details/useInventoryStockMoveItem'
import { useInventoryStockQuickAdjust } from '~/components/inventory/stock-details/useInventoryStockQuickAdjust'
import { useInventoryStockRecordUsage } from '~/components/inventory/stock-details/useInventoryStockRecordUsage'
import { INVENTORY_STOCKS_QUERY_KEY } from '~/types/inventory'
import { getErrorMessage } from '~/utils/errors'

type Props = {
  open: boolean
  stock: InventoryStockListItem | null
}

const props = defineProps<Props>()
const toast = useToast()
const queryClient = useQueryClient()
const inventoryStockStore = useInventoryStockStore()

const selectedMaterialId = computed<number>(() => {
  return props.stock?.material?.id ?? 0
})

const selectedMaterialQuery = useInventoryMaterialQuery(selectedMaterialId)

const isTogglingFavorite = computed<boolean>(() => {
  return inventoryStockStore.isMarkingFavorite || inventoryStockStore.isUnmarkingFavorite
})

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

const recordUsageProps = {
  get open() {
    return props.open
  },
  get stock() {
    return props.stock
  },
  get materialDetail() {
    return selectedMaterialQuery.data.value ?? null
  },
}

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
} = useInventoryStockRecordUsage(recordUsageProps)

const toggleFavorite = async (): Promise<void> => {
  const stock = props.stock
  if (!stock || isTogglingFavorite.value) return

  try {
    if (stock.is_favorite) {
      await inventoryStockStore.unmarkFavorite(stock.id)
      toast.add({
        title: 'Removed from favorites',
        color: 'success',
        duration: 2500,
      })
    } else {
      await inventoryStockStore.markFavorite(stock.id)
      toast.add({
        title: 'Marked as favorite',
        color: 'success',
        duration: 2500,
      })
    }

    await queryClient.invalidateQueries({ queryKey: INVENTORY_STOCKS_QUERY_KEY })
  } catch (err: unknown) {
    toast.add({
      title: 'Failed to update favorite state',
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  }
}

const archiveItem = async (): Promise<void> => {
  const stock = props.stock
  if (!stock || inventoryStockStore.isArchivingStock) return

  try {
    await inventoryStockStore.archiveStock(stock.id)
    await queryClient.invalidateQueries({ queryKey: INVENTORY_STOCKS_QUERY_KEY })
    toast.add({
      title: 'Stock item archived',
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: 'Failed to archive stock item',
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  }
}
</script>

<template>
  <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
    <InventoryStockQuickActionsBar
      :is-favorite="Boolean(props.stock?.is_favorite)"
      :is-editing-stock="isEditingStock"
      :is-moving-stock="isMovingStock"
      :is-recording-usage="isRecordingUsage"
      :is-toggling-favorite="isTogglingFavorite"
      :is-archiving-stock="inventoryStockStore.isArchivingStock"
      @adjust-stock="openEditMode"
      @move-item="openMoveMode"
      @record-usage="openUsageMode"
      @toggle-favorite="toggleFavorite"
      @archive-item="archiveItem"
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
