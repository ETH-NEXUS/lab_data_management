<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import { useInventoryMaterialsQuery } from '~/composables/inventory/useInventoryMaterialQuery'
import { useInventoryMaterialOrderUnitsQuery } from '~/composables/inventory/useInventoryMaterialUnitQuery'
import type { CreateInventoryOrderPayload, InventoryMaterialListItem } from '~/types/inventory'
import { getErrorMessage } from '~/utils/errors'

type Props = {
  open: boolean
  isSubmitting?: boolean
}

type OrderFormState = {
  materialId: string
  orderUnitId: string
  amount: string
  orderDate: string
  status: string
  notes: string
}

const props = withDefaults(defineProps<Props>(), { isSubmitting: false })

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
  (e: 'submit', payload: CreateInventoryOrderPayload): void
}>()

const toast = useToast()
const materialsQuery = useInventoryMaterialsQuery()

/**
 * Returns current local date as YYYY-MM-DD for `<input type="date">`.
 *
 * Returned data example:
 * - `'2026-04-08'`
 */
const getTodayDate = (): string => {
  const now = new Date()
  const year = String(now.getFullYear())
  const month = String(now.getMonth() + 1).padStart(2, '0')
  const day = String(now.getDate()).padStart(2, '0')
  return `${year}-${month}-${day}`
}

const buildInitialFormState = (): OrderFormState => ({
  materialId: '',
  orderUnitId: '',
  amount: '',
  orderDate: getTodayDate(),
  status: 'ordered',
  notes: '',
})

const formState = ref<OrderFormState>(buildInitialFormState())
const selectedMaterialIdRef = computed<number>(() => {
  const parsedId = Number.parseInt(formState.value.materialId, 10)
  return Number.isInteger(parsedId) && parsedId > 0 ? parsedId : 0
})
const materialOrderUnitsQuery = useInventoryMaterialOrderUnitsQuery(selectedMaterialIdRef)

const statusOptions = [
  { value: 'ordered', label: 'Ordered' },
  { value: 'tentative', label: 'Tentative' },
  { value: 'product_arrived', label: 'Product arrived' },
]

const orderUnitOptions = computed(() => {
  const materialUnits = materialOrderUnitsQuery.data.value ?? []
  return materialUnits.map((unit) => ({
    id: unit.id,
    label: unit.display_name || unit.unit?.label || unit.unit?.name || `Unit #${unit.id}`,
  }))
})

const isOrderUnitsLoading = computed<boolean>(() => {
  return selectedMaterialIdRef.value > 0 && materialOrderUnitsQuery.isPending.value
})

const orderUnitsErrorMessage = computed<string | null>(() => {
  const err = materialOrderUnitsQuery.error.value
  if (!err) {
    return null
  }
  return getErrorMessage(err)
})

const sortedMaterials = computed<InventoryMaterialListItem[]>(() => {
  const materials = [...(materialsQuery.data.value ?? [])]
  materials.sort((leftMaterial, rightMaterial) => {
    const leftLabel = (leftMaterial.label || leftMaterial.product_name || '').toLowerCase()
    const rightLabel = (rightMaterial.label || rightMaterial.product_name || '').toLowerCase()
    return leftLabel.localeCompare(rightLabel)
  })
  return materials
})

/**
 * Determines whether all required fields are valid for submission.
 *
 * Checked values:
 * - materialId: positive integer
 * - orderUnitId: positive integer
 * - amount: finite number > 0
 * - orderDate: non-empty YYYY-MM-DD
 * - status: non-empty
 */
const canSubmit = computed<boolean>(() => {
  if (props.isSubmitting || isOrderUnitsLoading.value) {
    return false
  }

  const materialId = Number.parseInt(formState.value.materialId, 10)
  const orderUnitId = Number.parseInt(formState.value.orderUnitId, 10)
  const amountValue = Number(formState.value.amount)
  const orderDate = formState.value.orderDate.trim()
  const status = formState.value.status.trim()

  return (
    Number.isInteger(materialId) &&
    materialId > 0 &&
    Number.isInteger(orderUnitId) &&
    orderUnitId > 0 &&
    Number.isFinite(amountValue) &&
    amountValue > 0 &&
    orderDate !== '' &&
    status !== ''
  )
})

watch(
  () => props.open,
  (isOpen) => {
    if (!isOpen) {
      return
    }

    formState.value = buildInitialFormState()
  },
)

watch(
  () => formState.value.materialId,
  (materialIdText) => {
    formState.value.orderUnitId = ''
    if (materialIdText.trim() === '') {
      return
    }
  },
)

const closeModal = (): void => {
  emit('update:open', false)
}

/**
 * Converts YYYY-MM-DD date input into ISO datetime accepted by order API.
 *
 * Input example:
 * - `'2026-04-08'`
 *
 * Returned value example:
 * - `'2026-04-08T00:00:00Z'`
 */
const toApiDateTime = (dateValue: string): string => `${dateValue}T00:00:00Z`

const submitForm = (): void => {
  const materialId = Number.parseInt(formState.value.materialId, 10)
  const orderUnitId = Number.parseInt(formState.value.orderUnitId, 10)
  const amountValue = Number(formState.value.amount)
  const orderDate = formState.value.orderDate.trim()
  const status = formState.value.status.trim()

  if (!Number.isInteger(materialId) || materialId <= 0) {
    toast.add({
      title: 'Material is required',
      color: 'warning',
      duration: 2500,
    })
    return
  }

  if (!Number.isInteger(orderUnitId) || orderUnitId <= 0) {
    toast.add({
      title: 'Order unit is required',
      color: 'warning',
      duration: 2500,
    })
    return
  }

  if (!Number.isFinite(amountValue) || amountValue <= 0) {
    toast.add({
      title: 'Amount must be greater than zero',
      color: 'warning',
      duration: 2500,
    })
    return
  }

  if (orderDate === '') {
    toast.add({
      title: 'Order date is required',
      color: 'warning',
      duration: 2500,
    })
    return
  }

  if (status === '') {
    toast.add({
      title: 'Status is required',
      color: 'warning',
      duration: 2500,
    })
    return
  }

  emit('submit', {
    material_id: materialId,
    order_unit_id: orderUnitId,
    amount: String(amountValue),
    order_date: toApiDateTime(orderDate),
    status,
    notes: formState.value.notes.trim() || null,
  })
}
</script>

<template>
  <UModal
    :open="props.open"
    title="Register order"
    description="Create a new inventory order entry."
    class="w-full sm:max-w-2xl"
    :ui="{ content: 'rounded-2xl bg-white shadow-md' }"
    @update:open="(isOpen) => emit('update:open', isOpen)"
  >
    <template #body>
      <div class="space-y-4 p-6">
        <div class="grid gap-4 sm:grid-cols-2">
          <div class="space-y-1">
            <label class="block text-sm font-medium text-slate-700">Material</label>
            <select
              v-model="formState.materialId"
              class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
              :disabled="materialsQuery.isPending.value || props.isSubmitting"
            >
              <option value="">
                {{ materialsQuery.isPending.value ? 'Loading materials...' : 'Select material' }}
              </option>
              <option v-for="material in sortedMaterials" :key="material.id" :value="String(material.id)">
                {{ material.label || material.product_name || `Material #${material.id}` }}
              </option>
            </select>
          </div>

          <div class="space-y-1">
            <label class="block text-sm font-medium text-slate-700">Order unit</label>
            <select
              v-model="formState.orderUnitId"
              class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50 disabled:cursor-not-allowed disabled:opacity-70"
              :disabled="isOrderUnitsLoading || formState.materialId === '' || props.isSubmitting"
            >
              <option value="">
                {{
                  isOrderUnitsLoading
                    ? 'Loading units...'
                    : orderUnitsErrorMessage
                      ? 'Failed to load units'
                      : orderUnitOptions.length > 0
                        ? 'Select order unit'
                        : 'No order units available'
                }}
              </option>
              <option v-for="unitOption in orderUnitOptions" :key="unitOption.id" :value="String(unitOption.id)">
                {{ unitOption.label }}
              </option>
            </select>
            <p v-if="orderUnitsErrorMessage" class="text-xs text-red-600">
              {{ orderUnitsErrorMessage }}
            </p>
          </div>

          <div class="space-y-1">
            <label class="block text-sm font-medium text-slate-700">Amount</label>
            <input
              v-model="formState.amount"
              type="number"
              min="0.000001"
              step="0.000001"
              class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
              :disabled="props.isSubmitting"
            />
          </div>

          <div class="space-y-1">
            <label class="block text-sm font-medium text-slate-700">Order date</label>
            <input
              v-model="formState.orderDate"
              type="date"
              class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
              :disabled="props.isSubmitting"
            />
          </div>

          <div class="space-y-1 sm:col-span-2">
            <label class="block text-sm font-medium text-slate-700">Status</label>
            <select
              v-model="formState.status"
              class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
              :disabled="props.isSubmitting"
            >
              <option v-for="statusOption in statusOptions" :key="statusOption.value" :value="statusOption.value">
                {{ statusOption.label }}
              </option>
            </select>
          </div>

          <div class="space-y-1 sm:col-span-2">
            <label class="block text-sm font-medium text-slate-700">Notes</label>
            <textarea
              v-model="formState.notes"
              rows="3"
              class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
              :disabled="props.isSubmitting"
            />
          </div>
        </div>
      </div>
    </template>

    <template #footer>
      <div class="flex w-full justify-end gap-2 px-6 pb-6">
        <UButton variant="ghost" color="neutral" label="Cancel" :disabled="props.isSubmitting" @click="closeModal" />
        <UButton
          color="primary"
          label="Register order"
          :loading="props.isSubmitting"
          :disabled="!canSubmit"
          @click="submitForm"
        />
      </div>
    </template>
  </UModal>
</template>
