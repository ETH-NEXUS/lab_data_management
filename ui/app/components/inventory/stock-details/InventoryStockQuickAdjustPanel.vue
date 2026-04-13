<script setup lang="ts">
import { computed } from 'vue'

type Props = {
  open: boolean
  isSavingStockAdjustment: boolean
  editQuantity: string
  editMinimumQuantity: string
  editNotes: string
  stockUnitLabel: string
}

const props = defineProps<Props>()

const emit = defineEmits<{
  (e: 'update:edit-quantity' | 'update:edit-minimum-quantity' | 'update:edit-notes', value: string): void
  (e: 'cancel' | 'save'): void
}>()

const { t } = useI18n()

const editQuantityModel = computed<string>({
  get: () => props.editQuantity,
  set: (value) => emit('update:edit-quantity', value),
})

const editMinimumQuantityModel = computed<string>({
  get: () => props.editMinimumQuantity,
  set: (value) => emit('update:edit-minimum-quantity', value),
})

const editNotesModel = computed<string>({
  get: () => props.editNotes,
  set: (value) => emit('update:edit-notes', value),
})
</script>

<template>
  <div v-if="props.open" class="space-y-3 rounded-lg border border-slate-200 bg-slate-50 p-3">
    <div class="grid gap-3 sm:grid-cols-2">
      <div class="space-y-1">
        <label class="block text-sm font-medium text-slate-700">{{
          t('inventory.stock_drawer.fields.quantity')
        }}</label>
        <input
          v-model="editQuantityModel"
          type="text"
          class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
          :disabled="props.isSavingStockAdjustment"
        />
      </div>
      <div class="space-y-1">
        <label class="block text-sm font-medium text-slate-700">{{
          t('inventory.stock_drawer.fields.minimum_quantity')
        }}</label>
        <input
          v-model="editMinimumQuantityModel"
          type="text"
          class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
          :disabled="props.isSavingStockAdjustment"
        />
      </div>
      <p class="text-xs text-slate-600 sm:col-span-2">
        {{ t('inventory.stock_drawer.fields.stock_unit') }}: {{ props.stockUnitLabel }}
      </p>
      <div class="space-y-1 sm:col-span-2">
        <label class="block text-sm font-medium text-slate-700">{{ t('inventory.stock_drawer.fields.notes') }}</label>
        <textarea
          v-model="editNotesModel"
          rows="3"
          class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
          :disabled="props.isSavingStockAdjustment"
        />
      </div>
    </div>
    <div class="flex justify-end gap-2">
      <UButton
        variant="ghost"
        color="neutral"
        :label="t('common.actions.cancel')"
        :disabled="props.isSavingStockAdjustment"
        @click="emit('cancel')"
      />
      <UButton
        color="primary"
        :label="t('common.actions.save')"
        :loading="props.isSavingStockAdjustment"
        :disabled="props.isSavingStockAdjustment"
        @click="emit('save')"
      />
    </div>
  </div>
</template>
