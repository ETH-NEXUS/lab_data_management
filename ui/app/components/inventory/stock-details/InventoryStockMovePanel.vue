<script setup lang="ts">
import { computed } from 'vue'
import type { InventoryRoom, InventorySector } from '~/types/inventory'

type Props = {
  open: boolean
  isSavingMove: boolean
  selectedRoomId: string
  selectedSectorId: string
  rooms: InventoryRoom[]
  filteredSectors: InventorySector[]
  isLookupsLoading: boolean
  lookupsErrorMessage?: string | null
  isMoveSaveDisabled: boolean
}

const props = withDefaults(defineProps<Props>(), {
  lookupsErrorMessage: null,
})

const emit = defineEmits<{
  (e: 'update:selected-room-id' | 'update:selected-sector-id', value: string): void
  (e: 'cancel' | 'save'): void
}>()

const { t } = useI18n()

const selectedRoomIdModel = computed<string>({
  get: () => props.selectedRoomId,
  set: (value) => emit('update:selected-room-id', value),
})

const selectedSectorIdModel = computed<string>({
  get: () => props.selectedSectorId,
  set: (value) => emit('update:selected-sector-id', value),
})
</script>

<template>
  <div v-if="props.open" class="space-y-3 rounded-lg border border-slate-200 bg-slate-50 p-3">
    <div class="grid gap-3 sm:grid-cols-2">
      <div class="space-y-1">
        <label class="block text-sm font-medium text-slate-700">{{ t('inventory.stock_drawer.fields.room') }}</label>
        <select
          v-model="selectedRoomIdModel"
          class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
          :disabled="props.isSavingMove || props.isLookupsLoading"
        >
          <option value="">
            {{
              props.isLookupsLoading
                ? t('inventory.stock_drawer.move.loading_rooms')
                : t('inventory.stock_drawer.move.select_room')
            }}
          </option>
          <option v-for="room in props.rooms" :key="room.id" :value="String(room.id)">
            {{ room.label || room.name }}
          </option>
        </select>
      </div>
      <div class="space-y-1">
        <label class="block text-sm font-medium text-slate-700">{{ t('inventory.stock_drawer.fields.sector') }}</label>
        <select
          v-model="selectedSectorIdModel"
          class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
          :disabled="props.isSavingMove || props.isLookupsLoading || props.selectedRoomId === ''"
        >
          <option value="">
            {{
              props.selectedRoomId === ''
                ? t('inventory.stock_drawer.move.select_room_first')
                : props.filteredSectors.length > 0
                  ? t('inventory.stock_drawer.move.select_sector')
                  : t('inventory.stock_drawer.move.no_sectors_available')
            }}
          </option>
          <option v-for="sector in props.filteredSectors" :key="sector.id" :value="String(sector.id)">
            {{ sector.label || sector.name }}
          </option>
        </select>
      </div>
    </div>

    <p v-if="props.lookupsErrorMessage" class="text-xs text-red-600">{{ props.lookupsErrorMessage }}</p>

    <div class="flex justify-end gap-2">
      <UButton
        variant="ghost"
        color="neutral"
        :label="t('common.actions.cancel')"
        :disabled="props.isSavingMove"
        @click="emit('cancel')"
      />
      <UButton
        color="primary"
        :label="t('common.actions.save')"
        :loading="props.isSavingMove"
        :disabled="props.isMoveSaveDisabled"
        @click="emit('save')"
      />
    </div>
  </div>
</template>
