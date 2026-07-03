<script setup lang="ts">
import { computed } from 'vue'
import type { InventoryRoom, InventorySector } from '~/types/inventory'

type Props = {
  open: boolean
  isSavingMove: boolean
  selectedRoomId: string
  selectedSectorIds: string[]
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
  (e: 'update:selected-room-id', value: string): void
  (e: 'update:selected-sector-ids', value: string[]): void
  (e: 'cancel' | 'save'): void
}>()

const { t } = useI18n()

const selectedRoomIdModel = computed<string>({
  get: () => props.selectedRoomId,
  set: (value) => emit('update:selected-room-id', value),
})

const selectedSectorIdsModel = computed<string[]>({
  get: () => props.selectedSectorIds,
  set: (value) => emit('update:selected-sector-ids', value),
})
</script>

<template>
  <div v-if="props.open" class="space-y-3 rounded-lg border border-slate-200 bg-slate-50 p-3">
    <div class="grid gap-3">
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
        <div
          class="max-h-40 space-y-2 overflow-y-auto rounded-md border border-slate-200 bg-white px-3 py-2"
          :class="{
            'cursor-not-allowed bg-slate-50 opacity-60': props.isSavingMove || props.isLookupsLoading || props.selectedRoomId === '',
          }"
        >
          <label
            v-for="sector in props.filteredSectors"
            :key="sector.id"
            class="flex items-center gap-2 text-sm text-slate-700"
          >
            <input
              v-model="selectedSectorIdsModel"
              type="checkbox"
              :value="String(sector.id)"
              :disabled="props.isSavingMove || props.isLookupsLoading || props.selectedRoomId === ''"
              class="h-4 w-4 rounded border-slate-300 text-blue-600 focus:ring-blue-500"
            >
            <span>{{ sector.label || sector.name }}</span>
          </label>
        </div>
        <p class="text-xs text-slate-500">
          {{
            props.selectedRoomId === ''
              ? t('inventory.stock_drawer.move.select_room_first')
              : props.filteredSectors.length > 0
                ? 'Select one or more sectors.'
                : t('inventory.stock_drawer.move.no_sectors_available')
          }}
        </p>
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
