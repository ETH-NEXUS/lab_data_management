<script setup lang="ts">
import { computed } from 'vue'
import InventoryAddItemValidationSection from '~/components/inventory/add-item/InventoryAddItemValidationSection.vue'
import { useInventoryAddItemForm } from '~/components/inventory/add-item/useInventoryAddItemForm'

type Props = {
  open: boolean
}

const props = defineProps<Props>()

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
}>()
const { t } = useI18n()

const closeModal = (): void => {
  emit('update:open', false)
}

const openRef = computed<boolean>(() => props.open)

const {
  formState,
  selectedOrderId,
  selectedItemTypeId,
  selectedMaterialId,
  selectedRoomId,
  roomOptions,
  supplierDetails,
  filteredMaterials,
  orderPrefillOptions,
  sectorOptions,
  stockUnitOptions,
  isStockUnitsLoading,
  materialsQuery,
  lookupsQuery,
  ordersQuery,
  materialsErrorMessage,
  sectorsErrorMessage,
  stockUnitsErrorMessage,
  ordersErrorMessage,
  validationMessages,
  isFormValid,
  isOrderPrefillDisabled,
  isSaveDisabled,
  isSubmitting,
  requestOrderPrefill,
  submitForm,
} = useInventoryAddItemForm({
  open: openRef,
  onSaved: closeModal,
})

/**
 * Builds sorted item-type select options from lookup data.
 *
 * Returned data example:
 * - `[{ id: 4, label: 'tube' }, { id: 5, label: 'plate' }]`
 */
const itemTypeOptions = computed<Array<{ id: number; label: string }>>(() => {
  const itemTypes = [...(lookupsQuery.data.value?.itemTypes ?? [])]

  itemTypes.sort((leftItemType, rightItemType) => {
    const leftLabel = (leftItemType.label || leftItemType.name || '').toLowerCase()
    const rightLabel = (rightItemType.label || rightItemType.name || '').toLowerCase()
    return leftLabel.localeCompare(rightLabel)
  })

  return itemTypes.map((itemType) => ({
    id: itemType.id,
    label: itemType.label || itemType.name || `Item type #${itemType.id}`,
  }))
})

const hasSelectedItemType = computed<boolean>(() => selectedItemTypeId.value > 0)
const hasSelectedMaterial = computed<boolean>(() => selectedMaterialId.value > 0)
const hasSelectedRoom = computed<boolean>(() => selectedRoomId.value > 0)
</script>

<template>
  <UModal
    :open="props.open"
    title="Add new item"
    description="Choose an item type first, optionally prefill from an order, then complete the stock fields."
    class="w-full sm:max-w-2xl"
    :ui="{ content: 'rounded-2xl bg-white shadow-md' }"
    @update:open="(isOpen) => emit('update:open', isOpen)"
  >
    <template #body>
      <div class="space-y-5 p-6">
        <section class="space-y-3">
          <p class="text-sm font-semibold text-slate-800">Item type</p>
          <div class="rounded-lg border border-slate-200 bg-slate-50 p-3">
            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">Item Type *</label>
              <select
                v-model="formState.itemTypeId"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                :disabled="lookupsQuery.isPending.value"
              >
                <option value="">
                  {{
                    lookupsQuery.isPending.value
                      ? 'Loading item types...'
                      : itemTypeOptions.length > 0
                        ? 'Select item type'
                        : 'No item types available'
                  }}
                </option>
                <option v-for="itemTypeOption in itemTypeOptions" :key="itemTypeOption.id" :value="String(itemTypeOption.id)">
                  {{ itemTypeOption.label }}
                </option>
              </select>
              <p class="text-xs text-slate-500">Step 1. Choose the item category to narrow materials and matching orders.</p>
            </div>
          </div>
        </section>

        <section class="space-y-3">
          <p class="text-sm font-semibold text-slate-800">{{ t('inventory.add_item.prefill.title') }}</p>
          <div class="rounded-lg border border-slate-200 bg-slate-50 p-3">
            <div class="grid gap-3 sm:grid-cols-[minmax(0,1fr)_auto] sm:items-end">
              <div class="space-y-1">
                <label class="block text-sm font-medium text-slate-700">
                  {{ t('inventory.add_item.prefill.order_label') }}
                </label>
                <select
                  v-model="selectedOrderId"
                  class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                  :disabled="ordersQuery.isPending.value || isSubmitting || !hasSelectedItemType"
                >
                  <option value="">
                    {{
                      !hasSelectedItemType
                        ? 'Select item type first'
                        : ordersQuery.isPending.value
                        ? t('inventory.add_item.prefill.select_options.loading_orders')
                        : ordersErrorMessage
                          ? t('inventory.add_item.prefill.select_options.failed_orders')
                          : orderPrefillOptions.length > 0
                            ? t('inventory.add_item.prefill.select_options.select_order')
                            : t('inventory.add_item.prefill.select_options.no_orders')
                    }}
                  </option>
                  <option
                    v-for="orderOption in orderPrefillOptions"
                    :key="orderOption.id"
                    :value="String(orderOption.id)"
                  >
                    {{ orderOption.label }}
                  </option>
                </select>
                <p v-if="ordersErrorMessage" class="text-xs text-red-600">
                  {{ ordersErrorMessage }}
                </p>
                <p v-else class="text-xs text-slate-500">Step 2. Optionally select an order from the chosen item type to populate fields.</p>
              </div>

              <UButton
                color="neutral"
                variant="soft"
                icon="i-heroicons-arrow-down-tray"
                :label="t('inventory.add_item.prefill.prefill_button')"
                :disabled="isOrderPrefillDisabled || !hasSelectedItemType"
                @click="requestOrderPrefill"
              />
            </div>
          </div>
        </section>

        <section class="space-y-3">
          <p class="text-sm font-semibold text-slate-800">Stock fields</p>
          <p class="text-xs text-slate-500">Step 3. Complete the remaining stock details for the selected material.</p>

          <div class="grid gap-4 sm:grid-cols-2">
            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">Material *</label>
              <select
                v-model="formState.materialId"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                :disabled="materialsQuery.isPending.value || !hasSelectedItemType"
              >
                <option value="">
                  {{
                    !hasSelectedItemType
                      ? 'Select item type first'
                      : materialsQuery.isPending.value
                      ? 'Loading materials...'
                      : materialsErrorMessage
                        ? 'Failed to load materials'
                        : filteredMaterials.length > 0
                          ? 'Select material'
                          : 'No materials available'
                  }}
                </option>
                <option v-for="material in filteredMaterials" :key="material.id" :value="String(material.id)">
                  {{ material.label || material.product_name || `Material #${material.id}` }}
                </option>
              </select>
              <p v-if="materialsErrorMessage" class="text-xs text-red-600">
                {{ materialsErrorMessage }}
              </p>
            </div>

            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">Room *</label>
              <select
                v-model="formState.roomId"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                :disabled="lookupsQuery.isPending.value"
              >
                <option value="">
                  {{
                    lookupsQuery.isPending.value
                      ? 'Loading rooms...'
                      : roomOptions.length > 0
                        ? 'Select room'
                        : 'No rooms available'
                  }}
                </option>
                <option v-for="roomOption in roomOptions" :key="roomOption.id" :value="String(roomOption.id)">
                  {{ roomOption.label }}
                </option>
              </select>
            </div>

            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">Sector *</label>
              <select
                v-model="formState.sectorId"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                :disabled="lookupsQuery.isPending.value || !hasSelectedRoom"
              >
                <option value="">
                  {{
                    !hasSelectedRoom
                      ? 'Select room first'
                      : lookupsQuery.isPending.value
                      ? 'Loading sectors...'
                      : sectorsErrorMessage
                        ? 'Failed to load sectors'
                        : sectorOptions.length > 0
                          ? 'Select sector'
                          : 'No sectors available'
                  }}
                </option>
                <option v-for="sectorOption in sectorOptions" :key="sectorOption.id" :value="String(sectorOption.id)">
                  {{ sectorOption.label }}
                </option>
              </select>
              <p v-if="sectorsErrorMessage" class="text-xs text-red-600">
                {{ sectorsErrorMessage }}
              </p>
            </div>

            <div class="space-y-2 sm:col-span-2">
              <p class="text-sm font-medium text-slate-700">Supplier and catalog</p>
              <div class="grid gap-3 rounded-lg border border-slate-200 bg-slate-50 p-3 sm:grid-cols-2">
                <div class="space-y-1">
                  <p class="text-xs font-medium uppercase tracking-wide text-slate-500">Manufacturer</p>
                  <p class="text-sm text-slate-700">
                    {{ hasSelectedMaterial ? (supplierDetails.manufacturer || '—') : 'Select material first' }}
                  </p>
                </div>

                <div class="space-y-1">
                  <p class="text-xs font-medium uppercase tracking-wide text-slate-500">Supplier</p>
                  <p class="text-sm text-slate-700">
                    {{ hasSelectedMaterial ? (supplierDetails.vendor || '—') : 'Select material first' }}
                  </p>
                </div>

                <div class="space-y-1">
                  <p class="text-xs font-medium uppercase tracking-wide text-slate-500">Manufacturer catalog number</p>
                  <p class="text-sm text-slate-700">
                    {{ hasSelectedMaterial ? (supplierDetails.manufacturerCatalogNumber || '—') : 'Select material first' }}
                  </p>
                </div>

                <div class="space-y-1">
                  <p class="text-xs font-medium uppercase tracking-wide text-slate-500">Supplier catalog number</p>
                  <p class="text-sm text-slate-700">
                    {{ hasSelectedMaterial ? (supplierDetails.vendorCatalogNumber || '—') : 'Select material first' }}
                  </p>
                </div>
              </div>
            </div>

            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">Stock unit *</label>
              <select
                v-model="formState.stockUnitId"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                :disabled="selectedMaterialId <= 0 || isStockUnitsLoading"
              >
                <option value="">
                  {{
                    selectedMaterialId <= 0
                      ? 'Select material first'
                      : isStockUnitsLoading
                        ? 'Loading stock units...'
                        : stockUnitsErrorMessage
                          ? 'Failed to load stock units'
                          : stockUnitOptions.length > 0
                            ? 'Select stock unit'
                            : 'No stock units available'
                  }}
                </option>
                <option
                  v-for="stockUnitOption in stockUnitOptions"
                  :key="stockUnitOption.id"
                  :value="String(stockUnitOption.id)"
                >
                  {{ stockUnitOption.label }}
                </option>
              </select>
              <p v-if="stockUnitsErrorMessage" class="text-xs text-red-600">
                {{ stockUnitsErrorMessage }}
              </p>
            </div>

            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">Quantity *</label>
              <input
                v-model="formState.quantity"
                type="number"
                min="0.000001"
                step="0.000001"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
              />
            </div>

            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">Minimum quantity *</label>
              <input
                v-model="formState.minimumQuantity"
                type="number"
                min="0"
                step="0.000001"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
              />
            </div>

            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">Lot number</label>
              <input
                v-model="formState.lotNumber"
                type="text"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
              />
            </div>

            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">Expiry date</label>
              <input
                v-model="formState.expiryDate"
                type="date"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
              />
            </div>

            <div class="flex items-end pb-1">
              <label class="inline-flex items-center gap-2 text-sm text-slate-700">
                <input v-model="formState.isFavorite" type="checkbox" />
                Mark as favorite
              </label>
            </div>

            <div class="space-y-1 sm:col-span-2">
              <label class="block text-sm font-medium text-slate-700">Notes</label>
              <textarea
                v-model="formState.notes"
                rows="3"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
              />
            </div>

            <InventoryAddItemValidationSection :validation-messages="validationMessages" />
          </div>
        </section>
      </div>
    </template>

    <template #footer>
      <div class="flex w-full justify-end gap-2 px-6 pb-6">
        <UButton variant="ghost" color="neutral" label="Cancel" :disabled="isSubmitting" @click="closeModal" />
        <UButton
          :label="isFormValid ? 'Save item' : 'Complete required fields'"
          color="primary"
          :loading="isSubmitting"
          :disabled="isSaveDisabled"
          @click="submitForm"
        />
      </div>
    </template>
  </UModal>
</template>
