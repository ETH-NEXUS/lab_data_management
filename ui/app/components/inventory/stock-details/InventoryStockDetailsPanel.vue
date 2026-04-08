<script setup lang="ts">
import type { InventoryStockDetailsPanelProps } from '~/components/inventory/stock-details/useInventoryStockDetailsPanel'
import { useInventoryStockDetailsPanel } from '~/components/inventory/stock-details/useInventoryStockDetailsPanel'
import { useInventoryStockQuickAdjust } from '~/components/inventory/stock-details/useInventoryStockQuickAdjust'
import { useInventoryStockMoveItem } from '~/components/inventory/stock-details/useInventoryStockMoveItem'
import { formatDateTime } from '~/utils/dateTime'

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

        <div v-if="isMovingStock" class="space-y-3 rounded-lg border border-slate-200 bg-slate-50 p-3">
          <div class="grid gap-3 sm:grid-cols-2">
            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">{{
                t('inventory.stock_drawer.fields.room')
              }}</label>
              <select
                v-model="selectedRoomId"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
                :disabled="isSavingMove || isLookupsLoading"
              >
                <option value="">
                  {{
                    isLookupsLoading
                      ? t('inventory.stock_drawer.move.loading_rooms')
                      : t('inventory.stock_drawer.move.select_room')
                  }}
                </option>
                <option v-for="room in rooms" :key="room.id" :value="String(room.id)">
                  {{ room.label || room.name }}
                </option>
              </select>
            </div>
            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">{{
                t('inventory.stock_drawer.fields.sector')
              }}</label>
              <select
                v-model="selectedSectorId"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
                :disabled="isSavingMove || isLookupsLoading || selectedRoomId === ''"
              >
                <option value="">
                  {{
                    selectedRoomId === ''
                      ? t('inventory.stock_drawer.move.select_room_first')
                      : filteredSectors.length > 0
                        ? t('inventory.stock_drawer.move.select_sector')
                        : t('inventory.stock_drawer.move.no_sectors_available')
                  }}
                </option>
                <option v-for="sector in filteredSectors" :key="sector.id" :value="String(sector.id)">
                  {{ sector.label || sector.name }}
                </option>
              </select>
            </div>
          </div>
          <p v-if="lookupsErrorMessage" class="text-xs text-red-600">{{ lookupsErrorMessage }}</p>
          <div class="flex justify-end gap-2">
            <UButton
              variant="ghost"
              color="neutral"
              :label="t('common.actions.cancel')"
              :disabled="isSavingMove"
              @click="cancelMoveMode"
            />
            <UButton
              color="primary"
              :label="t('common.actions.save')"
              :loading="isSavingMove"
              :disabled="isMoveSaveDisabled"
              @click="saveMove"
            />
          </div>
        </div>
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

      <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
        <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
          {{ t('inventory.stock_drawer.sections.usage') }}
        </p>

        <div class="space-y-2 rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
          <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">
            {{ t('inventory.stock_drawer.fields.recent_usages') }}
          </p>

          <p v-if="props.isUsagesLoading" class="text-sm text-slate-500">
            {{ t('inventory.stock_workspace.loading') }}
          </p>
          <p v-else-if="sortedUsageEntries.length === 0" class="text-sm text-slate-600">
            {{ t('inventory.stock_drawer.values.no_usage') }}
          </p>
          <ul v-else class="space-y-1.5 text-sm text-slate-700">
            <li v-for="usage in sortedUsageEntries.slice(0, 8)" :key="usage.id">
              <span class="font-medium">
                {{
                  formatDateTime(
                    usage.used_at,
                    { dateStyle: 'medium', timeStyle: 'short' },
                    t('inventory.stock_drawer.values.none'),
                  )
                }}
              </span>
              · {{ usage.quantity_used }} {{ usage.usage_unit.display_name }}
              <span class="text-slate-500">
                ({{ usage.project?.label || usage.project?.name || t('inventory.stock_drawer.values.none') }})
              </span>
            </li>
          </ul>
        </div>

        <div class="space-y-2 rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
          <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">
            {{ t('inventory.stock_drawer.fields.linked_experiments_projects') }}
          </p>
          <p v-if="linkedExperimentProjectLabels.length === 0" class="text-sm text-slate-600">
            {{ t('inventory.stock_drawer.values.no_linked_items') }}
          </p>
          <ul v-else class="list-disc space-y-1 pl-5 text-sm text-slate-700">
            <li v-for="linkedLabel in linkedExperimentProjectLabels" :key="linkedLabel">
              {{ linkedLabel }}
            </li>
          </ul>
        </div>
      </section>

      <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
        <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
          {{ t('inventory.stock_drawer.sections.order') }}
        </p>

        <div class="space-y-2 rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
          <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">
            {{ t('inventory.stock_drawer.fields.order_history') }}
          </p>

          <p v-if="props.isOrdersLoading" class="text-sm text-slate-500">
            {{ t('inventory.stock_workspace.loading') }}
          </p>
          <p v-else-if="sortedOrderEntries.length === 0" class="text-sm text-slate-600">
            {{ t('inventory.stock_drawer.values.no_order_history') }}
          </p>
          <ul v-else class="space-y-1.5 text-sm text-slate-700">
            <li v-for="order in sortedOrderEntries.slice(0, 8)" :key="order.id">
              <span class="font-medium">
                {{ formatDateTime(order.order_date, { dateStyle: 'medium' }, t('inventory.stock_drawer.values.none')) }}
              </span>
              · {{ order.amount }} {{ order.order_unit.display_name }}
              <span class="text-slate-500"> ({{ order.status_label || order.status }}) </span>
            </li>
          </ul>
        </div>

        <div class="space-y-2 rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
          <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">
            {{ t('inventory.stock_drawer.fields.responsible_user') }}
          </p>
          <p v-if="responsibleUsers.length === 0" class="text-sm text-slate-600">
            {{ t('inventory.stock_drawer.values.no_responsible_user') }}
          </p>
          <ul v-else class="list-disc space-y-1 pl-5 text-sm text-slate-700">
            <li v-for="userLabel in responsibleUsers" :key="userLabel">
              {{ userLabel }}
            </li>
          </ul>
        </div>
      </section>

      <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
        <div class="flex items-center justify-between gap-3">
          <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
            {{ t('inventory.stock_drawer.sections.metadata') }}
          </p>

          <UButton
            size="xs"
            variant="ghost"
            color="neutral"
            :label="
              isMetadataExpanded
                ? t('inventory.stock_drawer.metadata_toggle.hide')
                : t('inventory.stock_drawer.metadata_toggle.show')
            "
            @click="toggleMetadata"
          />
        </div>

        <div v-if="isMetadataExpanded" class="grid gap-2 sm:grid-cols-2">
          <div
            v-for="field in metadataFields"
            :key="field.label"
            class="rounded-lg border border-slate-200 bg-slate-50 px-3 py-2"
          >
            <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">{{ field.label }}</p>
            <p class="mt-1 text-sm text-slate-800">{{ field.value }}</p>
          </div>
        </div>
      </section>
    </div>
  </aside>
</template>
