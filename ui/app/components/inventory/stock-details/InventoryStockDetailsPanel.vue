<script setup lang="ts">
import type { InventoryStockDetailsPanelProps } from '~/components/inventory/stock-details/useInventoryStockDetailsPanel'
import InventoryStockMovePanel from '~/components/inventory/stock-details/InventoryStockMovePanel.vue'
import InventoryStockMetadataSection from '~/components/inventory/stock-details/InventoryStockMetadataSection.vue'
import InventoryStockOrderSection from '~/components/inventory/stock-details/InventoryStockOrderSection.vue'
import InventoryStockRecordUsagePanel from '~/components/inventory/stock-details/InventoryStockRecordUsagePanel.vue'
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
    <header class="flex items-start justify-between gap-4 border-b border-[var(--app-border)] px-5 py-4">
      <div class="min-w-0 space-y-1">
        <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
          {{ t('inventory.stock_drawer.title') }}
        </p>
        <h2 class="truncate text-lg font-semibold text-slate-900">
          {{ props.stock.material.product_name }}
        </h2>
        <p class="text-xs text-slate-600">
          {{ t('inventory.stock_drawer.subtitle', { id: props.stock.id }) }}
        </p>
      </div>

      <UButton
        size="xs"
        variant="ghost"
        color="neutral"
        icon="i-heroicons-x-mark"
        :label="t('inventory.stock_drawer.actions.close')"
        @click="emit('close')"
      />
    </header>

    <div class="space-y-5 overflow-y-auto px-5 py-4">
      <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
        <div class="flex items-center justify-between gap-3">
          <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
            {{ t('inventory.page.section_labels.quick_actions') }}
          </p>
          <div class="flex items-center gap-2">
            <UButton
              color="primary"
              size="xs"
              :label="t('inventory.stock_drawer.actions.adjust_stock')"
              :disabled="isEditingStock"
              @click="openEditMode"
            />
            <UButton
              color="neutral"
              variant="soft"
              size="xs"
              :label="t('inventory.stock_drawer.actions.move_item')"
              :disabled="isMovingStock"
              @click="openMoveMode"
            />
            <UButton
              color="neutral"
              variant="soft"
              size="xs"
              label="Record usage"
              :disabled="isRecordingUsage"
              @click="openUsageMode"
            />
          </div>
        </div>

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

      <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
        <div class="flex flex-wrap items-center gap-2">
          <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
            {{ t('inventory.stock_drawer.sections.operational') }}
          </p>
          <UBadge :color="inventoryStatusColor" variant="soft">{{ inventoryStatusLabel }}</UBadge>
        </div>

        <div class="grid gap-2 sm:grid-cols-2">
          <div
            v-for="field in operationalFields"
            :key="field.label"
            :class="['rounded-lg border border-slate-200 bg-slate-50 px-3 py-2', field.wide ? 'sm:col-span-2' : '']"
          >
            <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">{{ field.label }}</p>
            <p class="mt-1 text-sm text-slate-800">{{ field.value }}</p>
          </div>
        </div>
      </section>

      <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
        <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
          {{ t('inventory.stock_drawer.sections.unit_conversion') }}
        </p>

        <div class="grid gap-2 sm:grid-cols-3">
          <div
            v-for="field in unitConversionFields"
            :key="field.label"
            class="rounded-lg border border-slate-200 bg-slate-50 px-3 py-2"
          >
            <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">{{ field.label }}</p>
            <p class="mt-1 text-sm text-slate-800">{{ field.value }}</p>
          </div>
        </div>
      </section>

      <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
        <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
          {{ t('inventory.stock_drawer.sections.material_identity') }}
        </p>

        <p v-if="props.isMaterialLoading" class="text-sm text-slate-500">
          {{ t('inventory.stock_workspace.loading') }}
        </p>

        <div v-else class="grid gap-2 sm:grid-cols-2">
          <div
            v-for="field in materialIdentityFields"
            :key="field.label"
            :class="['rounded-lg border border-slate-200 bg-slate-50 px-3 py-2', field.wide ? 'sm:col-span-2' : '']"
          >
            <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">{{ field.label }}</p>
            <p class="mt-1 text-sm text-slate-800">{{ field.value }}</p>
          </div>
        </div>
      </section>

      <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
        <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
          {{ t('inventory.stock_drawer.sections.supplier') }}
        </p>

        <div class="grid gap-2 sm:grid-cols-2">
          <div
            v-for="field in supplierFields"
            :key="field.label"
            :class="['rounded-lg border border-slate-200 bg-slate-50 px-3 py-2', field.wide ? 'sm:col-span-2' : '']"
          >
            <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">{{ field.label }}</p>
            <p class="mt-1 text-sm text-slate-800">{{ field.value }}</p>
          </div>
        </div>
      </section>

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
