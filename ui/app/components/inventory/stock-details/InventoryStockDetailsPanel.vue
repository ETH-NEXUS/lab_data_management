<script setup lang="ts">
import type { InventoryStockDetailsPanelProps } from '~/components/inventory/stock-details/useInventoryStockDetailsPanel'
import InventoryStockDetailsHeader from '~/components/inventory/stock-details/InventoryStockDetailsHeader.vue'
import InventoryStockMaterialIdentitySection from '~/components/inventory/stock-details/InventoryStockMaterialIdentitySection.vue'
import InventoryStockMovePanel from '~/components/inventory/stock-details/InventoryStockMovePanel.vue'
import InventoryStockMetadataSection from '~/components/inventory/stock-details/InventoryStockMetadataSection.vue'
import InventoryStockOperationalSection from '~/components/inventory/stock-details/InventoryStockOperationalSection.vue'
import InventoryStockOrderSection from '~/components/inventory/stock-details/InventoryStockOrderSection.vue'
import InventoryStockQuickActionsBar from '~/components/inventory/stock-details/InventoryStockQuickActionsBar.vue'
import InventoryStockRecordUsagePanel from '~/components/inventory/stock-details/InventoryStockRecordUsagePanel.vue'
import InventoryStockSupplierSection from '~/components/inventory/stock-details/InventoryStockSupplierSection.vue'
import InventoryStockUnitConversionSection from '~/components/inventory/stock-details/InventoryStockUnitConversionSection.vue'
import InventoryStockUsageSection from '~/components/inventory/stock-details/InventoryStockUsageSection.vue'
import { useInventoryStockDetailsPanel } from '~/components/inventory/stock-details/useInventoryStockDetailsPanel'
import { useInventoryStockQuickAdjust } from '~/components/inventory/stock-details/useInventoryStockQuickAdjust'
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

const { t } = useI18n()

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
      <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
        <InventoryStockQuickActionsBar
          :is-editing-stock="isEditingStock"
          :is-moving-stock="isMovingStock"
          :is-recording-usage="isRecordingUsage"
          @adjust-stock="openEditMode"
          @move-item="openMoveMode"
          @record-usage="openUsageMode"
        />

        <div v-if="isEditingStock" class="space-y-3 rounded-lg border border-slate-200 bg-slate-50 p-3">
          <div class="grid gap-3 sm:grid-cols-2">
            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">{{
                t('inventory.stock_drawer.fields.quantity')
              }}</label>
              <input
                v-model="editQuantity"
                type="text"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
                :disabled="isSavingStockAdjustment"
              />
            </div>
            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">{{
                t('inventory.stock_drawer.fields.minimum_quantity')
              }}</label>
              <input
                v-model="editMinimumQuantity"
                type="text"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
                :disabled="isSavingStockAdjustment"
              />
            </div>
            <p class="text-xs text-slate-600 sm:col-span-2">
              {{ t('inventory.stock_drawer.fields.stock_unit') }}: {{ stockUnitLabel }}
            </p>
            <div class="space-y-1 sm:col-span-2">
              <label class="block text-sm font-medium text-slate-700">{{
                t('inventory.stock_drawer.fields.notes')
              }}</label>
              <textarea
                v-model="editNotes"
                rows="3"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
                :disabled="isSavingStockAdjustment"
              />
            </div>
          </div>
          <div class="flex justify-end gap-2">
            <UButton
              variant="ghost"
              color="neutral"
              :label="t('common.actions.cancel')"
              :disabled="isSavingStockAdjustment"
              @click="cancelEditMode"
            />
            <UButton
              color="primary"
              :label="t('common.actions.save')"
              :loading="isSavingStockAdjustment"
              :disabled="isSavingStockAdjustment"
              @click="saveStockAdjustment"
            />
          </div>
        </div>

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
