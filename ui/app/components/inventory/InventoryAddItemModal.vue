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
  filteredMaterials,
  orderPrefillOptions,
  sectorOptions,
  brandOptions,
  manufacturerOptions,
  vendorOptions,
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
  isSelectedMaterialReagent,
  requestOrderPrefill,
  setReagentSafetyDataSheet,
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

const onSafetyDataSheetChange = (event: Event): void => {
  const inputElement = event.target as HTMLInputElement | null
  const selectedFile = inputElement?.files?.[0] ?? null
  setReagentSafetyDataSheet(selectedFile)
}
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
                <option
                  v-for="itemTypeOption in itemTypeOptions"
                  :key="itemTypeOption.id"
                  :value="String(itemTypeOption.id)"
                >
                  {{ itemTypeOption.label }}
                </option>
              </select>
              <p class="text-xs text-slate-500">
                Step 1. Choose the item category to narrow materials and matching orders.
              </p>
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
                <p v-else class="text-xs text-slate-500">
                  Step 2. Optionally select an order from the chosen item type to populate fields.
                </p>
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
              <div
                class="max-h-40 space-y-2 overflow-y-auto rounded-md border border-slate-200 bg-white px-3 py-2"
                :class="{
                  'cursor-not-allowed bg-slate-50 opacity-60': lookupsQuery.isPending.value || !hasSelectedRoom,
                }"
              >
                <label
                  v-for="sectorOption in sectorOptions"
                  :key="sectorOption.id"
                  class="flex items-center gap-2 text-sm text-slate-700"
                >
                  <input
                    v-model="formState.sectorIds"
                    type="checkbox"
                    :value="String(sectorOption.id)"
                    :disabled="lookupsQuery.isPending.value || !hasSelectedRoom"
                    class="h-4 w-4 rounded border-slate-300 text-blue-600 focus:ring-blue-500"
                  />
                  <span>{{ sectorOption.label }}</span>
                </label>
              </div>
              <p class="text-xs text-slate-500">
                {{
                  !hasSelectedRoom
                    ? 'Select room first'
                    : lookupsQuery.isPending.value
                      ? 'Loading sectors...'
                      : sectorsErrorMessage
                        ? 'Failed to load sectors'
                        : sectorOptions.length > 0
                          ? 'Select one or more sectors.'
                          : 'No sectors available'
                }}
              </p>
              <p v-if="sectorsErrorMessage" class="text-xs text-red-600">
                {{ sectorsErrorMessage }}
              </p>
            </div>

            <div v-if="hasSelectedMaterial" class="space-y-2 sm:col-span-2">
              <p class="text-sm font-medium text-slate-700">Additional material details</p>
              <p class="text-xs text-slate-500">
                Optional. Fill in anything missing on this material so you don't have to edit it later from the table.
              </p>
              <div class="grid gap-3 rounded-lg border border-slate-200 bg-slate-50 p-3 sm:grid-cols-2">
                <div class="space-y-1">
                  <label class="block text-sm font-medium text-slate-700">Brand</label>
                  <select
                    v-model="formState.additionalBrandId"
                    class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                    :disabled="lookupsQuery.isPending.value"
                  >
                    <option value="">No brand</option>
                    <option v-for="brandOption in brandOptions" :key="brandOption.id" :value="String(brandOption.id)">
                      {{ brandOption.label }}
                    </option>
                  </select>
                </div>

                <div class="space-y-1">
                  <label class="block text-sm font-medium text-slate-700">Default cost</label>
                  <input
                    v-model="formState.additionalDefaultCost"
                    type="text"
                    inputmode="decimal"
                    placeholder="e.g. 12.50"
                    class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                  />
                </div>

                <div class="space-y-1">
                  <label class="block text-sm font-medium text-slate-700">Manufacturer</label>
                  <select
                    v-model="formState.additionalManufacturerId"
                    class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                    :disabled="lookupsQuery.isPending.value"
                  >
                    <option value="">No manufacturer</option>
                    <option
                      v-for="manufacturerOption in manufacturerOptions"
                      :key="manufacturerOption.id"
                      :value="String(manufacturerOption.id)"
                    >
                      {{ manufacturerOption.label }}
                    </option>
                  </select>
                </div>

                <div class="space-y-1">
                  <label class="block text-sm font-medium text-slate-700">Supplier</label>
                  <select
                    v-model="formState.additionalVendorId"
                    class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                    :disabled="lookupsQuery.isPending.value"
                  >
                    <option value="">No supplier</option>
                    <option
                      v-for="vendorOption in vendorOptions"
                      :key="vendorOption.id"
                      :value="String(vendorOption.id)"
                    >
                      {{ vendorOption.label }}
                    </option>
                  </select>
                </div>

                <div class="space-y-1">
                  <label class="block text-sm font-medium text-slate-700">Manufacturer catalog number</label>
                  <input
                    v-model="formState.additionalManufacturerCatalogNumber"
                    type="text"
                    class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                  />
                </div>

                <div class="space-y-1">
                  <label class="block text-sm font-medium text-slate-700">Supplier catalog number</label>
                  <input
                    v-model="formState.additionalVendorCatalogNumber"
                    type="text"
                    class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                  />
                </div>

                <div class="space-y-1">
                  <label class="block text-sm font-medium text-slate-700">Capacity value</label>
                  <input
                    v-model="formState.additionalCapacityValue"
                    type="text"
                    inputmode="decimal"
                    placeholder="e.g. 200"
                    class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                  />
                </div>

                <div class="space-y-1">
                  <label class="block text-sm font-medium text-slate-700">Capacity unit</label>
                  <input
                    v-model="formState.additionalCapacityUnit"
                    type="text"
                    placeholder="e.g. ul, ml, mm"
                    class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                  />
                </div>

                <div class="space-y-1 sm:col-span-2">
                  <label class="block text-sm font-medium text-slate-700">Description</label>
                  <textarea
                    v-model="formState.additionalDescription"
                    rows="2"
                    class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                  />
                </div>

                <div class="space-y-1">
                  <label class="block text-sm font-medium text-slate-700">Serial number</label>
                  <input
                    v-model="formState.additionalSerialNumber"
                    type="text"
                    class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                  />
                </div>

                <div class="space-y-1">
                  <label class="block text-sm font-medium text-slate-700">Order number</label>
                  <input
                    v-model="formState.additionalOrderNumber"
                    type="text"
                    class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                  />
                </div>

                <div class="space-y-1">
                  <label class="block text-sm font-medium text-slate-700">Shelf life (days)</label>
                  <input
                    v-model="formState.additionalLifetimeDays"
                    type="number"
                    min="0"
                    step="1"
                    inputmode="numeric"
                    class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                  />
                </div>

                <div class="space-y-1">
                  <label class="block text-sm font-medium text-slate-700">Active</label>
                  <select
                    v-model="formState.additionalIsActive"
                    class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                  >
                    <option value="">Don't change</option>
                    <option value="true">Active</option>
                    <option value="false">Inactive</option>
                  </select>
                </div>
              </div>
            </div>

            <div v-if="isSelectedMaterialReagent" class="space-y-2 sm:col-span-2">
              <p class="text-sm font-medium text-slate-700">Reagent metadata</p>
              <div class="grid gap-3 rounded-lg border border-slate-200 bg-slate-50 p-3 sm:grid-cols-2">
                <div class="space-y-1">
                  <label class="block text-sm font-medium text-slate-700">Storage temperature *</label>
                  <select
                    v-model="formState.reagentStorageTemperature"
                    class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                  >
                    <option value="">Select storage temperature</option>
                    <option value="4°C">4°C</option>
                    <option value="RT">RT</option>
                    <option value="-20°C">-20°C</option>
                    <option value="-80°C">-80°C</option>
                    <option value="LN">LN</option>
                  </select>
                </div>

                <div class="space-y-1">
                  <label class="block text-sm font-medium text-slate-700">Safety data sheet</label>
                  <input
                    type="file"
                    accept=".pdf,.doc,.docx"
                    class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                    @change="onSafetyDataSheetChange"
                  />
                  <p class="text-xs text-slate-500">
                    {{
                      formState.reagentSafetyDataSheet
                        ? `Selected file: ${formState.reagentSafetyDataSheet.name}`
                        : 'Optional. Upload one safety data sheet file.'
                    }}
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
                type="text"
                inputmode="decimal"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
              />
            </div>

            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">Minimum quantity</label>
              <input
                v-model="formState.minimumQuantity"
                type="number"
                min="0"
                step="1"
                placeholder="0"
                inputmode="numeric"
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
