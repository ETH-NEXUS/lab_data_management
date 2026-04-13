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

const closeModal = (): void => {
  emit('update:open', false)
}

const openRef = computed<boolean>(() => props.open)

const {
  formState,
  selectedOrderId,
  selectedMaterialId,
  sortedMaterials,
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
</script>

<template>
  <UModal
    :open="props.open"
    title="Add new item"
    description="Create one new inventory stock row using model-aligned fields."
    class="w-full sm:max-w-2xl"
    :ui="{ content: 'rounded-2xl bg-white shadow-md' }"
    @update:open="(isOpen) => emit('update:open', isOpen)"
  >
    <template #body>
      <div class="space-y-5 p-6">
        <section class="space-y-3">
          <p class="text-sm font-semibold text-slate-800">Prefill from order</p>
          <div class="rounded-lg border border-slate-200 bg-slate-50 p-3">
            <div class="grid gap-3 sm:grid-cols-[minmax(0,1fr)_auto] sm:items-end">
              <div class="space-y-1">
                <label class="block text-sm font-medium text-slate-700">Order</label>
                <select
                  v-model="selectedOrderId"
                  class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                  :disabled="ordersQuery.isPending.value || isSubmitting"
                >
                  <option value="">
                    {{
                      ordersQuery.isPending.value
                        ? 'Loading orders...'
                        : ordersErrorMessage
                          ? 'Failed to load orders'
                          : orderPrefillOptions.length > 0
                            ? 'Select order'
                            : 'No orders available'
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
              </div>

              <UButton
                color="neutral"
                variant="soft"
                icon="i-heroicons-arrow-down-tray"
                label="Prefill from order"
                :disabled="isOrderPrefillDisabled"
                @click="requestOrderPrefill"
              />
            </div>
          </div>
        </section>

        <section class="space-y-3">
          <p class="text-sm font-semibold text-slate-800">Stock fields</p>

          <div class="grid gap-4 sm:grid-cols-2">
            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">Material *</label>
              <select
                v-model="formState.materialId"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                :disabled="materialsQuery.isPending.value"
              >
                <option value="">
                  {{
                    materialsQuery.isPending.value
                      ? 'Loading materials...'
                      : materialsErrorMessage
                        ? 'Failed to load materials'
                        : sortedMaterials.length > 0
                          ? 'Select material'
                          : 'No materials available'
                  }}
                </option>
                <option v-for="material in sortedMaterials" :key="material.id" :value="String(material.id)">
                  {{ material.label || material.product_name || `Material #${material.id}` }}
                </option>
              </select>
              <p v-if="materialsErrorMessage" class="text-xs text-red-600">
                {{ materialsErrorMessage }}
              </p>
            </div>

            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">Sector *</label>
              <select
                v-model="formState.sectorId"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
                :disabled="lookupsQuery.isPending.value"
              >
                <option value="">
                  {{
                    lookupsQuery.isPending.value
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
