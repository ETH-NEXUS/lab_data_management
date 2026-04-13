<script setup lang="ts">
import { ref, watch } from 'vue'

type Props = {
  open: boolean
}

type AddItemFlow = 'existing_material_stock' | 'new_material_and_stock'

type AddItemFormState = {
  materialId: string
  sectorId: string
  stockUnitId: string
  quantity: string
  minimumQuantity: string
  lotNumber: string
  expiryDate: string
  notes: string
  isFavorite: boolean
}

const props = defineProps<Props>()

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
}>()

/**
 * Builds empty draft values for the add-item stock form.
 *
 * Returned data example:
 * - `{ materialId: '', sectorId: '', stockUnitId: '', quantity: '', minimumQuantity: '', lotNumber: '', expiryDate: '', notes: '', isFavorite: false }`
 */
const buildInitialFormState = (): AddItemFormState => ({
  materialId: '',
  sectorId: '',
  stockUnitId: '',
  quantity: '',
  minimumQuantity: '',
  lotNumber: '',
  expiryDate: '',
  notes: '',
  isFavorite: false,
})

const selectedFlow = ref<AddItemFlow>('existing_material_stock')
const formState = ref<AddItemFormState>(buildInitialFormState())

watch(
  () => props.open,
  (isOpen) => {
    if (!isOpen) {
      return
    }

    selectedFlow.value = 'existing_material_stock'
    formState.value = buildInitialFormState()
  },
)

const closeModal = (): void => {
  emit('update:open', false)
}
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
        <section class="space-y-2">
          <p class="text-sm font-semibold text-slate-800">Flow options</p>

          <div class="grid gap-2 sm:grid-cols-2">
            <button
              type="button"
              :class="[
                'rounded-lg border px-4 py-3 text-left',
                selectedFlow === 'existing_material_stock'
                  ? 'border-blue-300 bg-blue-50'
                  : 'border-slate-200 bg-slate-50 hover:border-slate-300',
              ]"
              @click="selectedFlow = 'existing_material_stock'"
            >
              <p class="text-sm font-semibold text-slate-800">Add stock for existing material</p>
              <p class="mt-1 text-sm text-slate-600">Creates one new row in the items table.</p>
            </button>

            <button
              type="button"
              class="cursor-not-allowed rounded-lg border border-slate-200 bg-slate-100 px-4 py-3 text-left opacity-70"
              disabled
            >
              <p class="text-sm font-semibold text-slate-700">Create new material and stock</p>
              <p class="mt-1 text-sm text-slate-600">Will be enabled after the stock-entry flow is finished.</p>
            </button>
          </div>
        </section>

        <section v-if="selectedFlow === 'existing_material_stock'" class="space-y-3">
          <p class="text-sm font-semibold text-slate-800">Stock fields</p>

          <div class="grid gap-4 sm:grid-cols-2">
            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">Material *</label>
              <select
                v-model="formState.materialId"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
              >
                <option value="">Select material (wired in next step)</option>
              </select>
            </div>

            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">Sector *</label>
              <select
                v-model="formState.sectorId"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
              >
                <option value="">Select sector (wired in next step)</option>
              </select>
            </div>

            <div class="space-y-1">
              <label class="block text-sm font-medium text-slate-700">Stock unit *</label>
              <select
                v-model="formState.stockUnitId"
                class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700"
              >
                <option value="">Select stock unit (wired in next step)</option>
              </select>
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
          </div>
        </section>
      </div>
    </template>

    <template #footer>
      <div class="flex w-full justify-end gap-2 px-6 pb-6">
        <UButton variant="ghost" color="neutral" label="Cancel" @click="closeModal" />
        <UButton color="primary" label="Save item" disabled />
      </div>
    </template>
  </UModal>
</template>
