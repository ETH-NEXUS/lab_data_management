<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import type { InventoryDashboardTilePreference } from '~/types/inventory'

type Props = {
  tiles: InventoryDashboardTilePreference[]
  isLoading: boolean
  isSaving: boolean
}

const props = defineProps<Props>()

const emit = defineEmits<{
  (e: 'save', tileKeys: string[]): void
}>()

const { t } = useI18n()
const isOpen = ref(false)
const selectedTileKeys = ref<string[]>([])
const minimumTileCount = 4
const maximumTileCount = 6

const sortedTiles = computed<InventoryDashboardTilePreference[]>(() => {
  return [...props.tiles].sort((leftTile, rightTile) => leftTile.position - rightTile.position)
})

const selectedTileCount = computed<number>(() => selectedTileKeys.value.length)
const canSave = computed<boolean>(() => {
  return selectedTileCount.value >= minimumTileCount && selectedTileCount.value <= maximumTileCount
})

watch(
  () => props.tiles,
  (tiles) => {
    selectedTileKeys.value = tiles
      .filter((tile) => tile.is_visible)
      .sort((leftTile, rightTile) => leftTile.position - rightTile.position)
      .map((tile) => tile.key)
  },
  { immediate: true },
)

const isTileSelected = (tileKey: string): boolean => selectedTileKeys.value.includes(tileKey)

const isTileDisabled = (tileKey: string): boolean => {
  if (isTileSelected(tileKey)) {
    return selectedTileCount.value <= minimumTileCount
  }

  return selectedTileCount.value >= maximumTileCount
}

const toggleTile = (tileKey: string, isSelected: boolean): void => {
  if (isSelected && !isTileSelected(tileKey) && selectedTileCount.value < maximumTileCount) {
    selectedTileKeys.value.push(tileKey)
  }

  if (!isSelected && isTileSelected(tileKey) && selectedTileCount.value > minimumTileCount) {
    selectedTileKeys.value = selectedTileKeys.value.filter((selectedTileKey) => selectedTileKey !== tileKey)
  }
}

const saveSelection = (): void => {
  if (!canSave.value) return
  emit('save', selectedTileKeys.value)
}
</script>

<template>
  <div class="flex justify-end">
    <UButton
      variant="ghost"
      color="neutral"
      icon="i-heroicons-adjustments-horizontal"
      :label="t('inventory.dashboard_tiles.configure')"
      :disabled="props.isLoading"
      @click="isOpen = true"
    />

    <UModal
      :open="isOpen"
      :title="t('inventory.dashboard_tiles.title')"
      :description="t('inventory.dashboard_tiles.description')"
      class="w-full sm:max-w-xl"
      :ui="{ content: 'rounded-2xl bg-white shadow-md' }"
      @update:open="(isModalOpen) => (isOpen = isModalOpen)"
    >
      <template #body>
        <div class="space-y-4 p-6">
          <p class="text-sm text-slate-600">
            {{ t('inventory.dashboard_tiles.selected_count', { count: selectedTileCount, maximum: maximumTileCount }) }}
          </p>

          <div class="space-y-2">
            <label
              v-for="tile in sortedTiles"
              :key="tile.key"
              class="flex cursor-pointer items-center gap-3 rounded-lg border border-slate-200 px-3 py-2 hover:bg-slate-50"
              :class="{ 'cursor-not-allowed opacity-60': isTileDisabled(tile.key) }"
            >
              <UCheckbox
                :model-value="isTileSelected(tile.key)"
                :disabled="props.isSaving || isTileDisabled(tile.key)"
                @update:model-value="(isSelected) => toggleTile(tile.key, Boolean(isSelected))"
              />
              <span class="text-sm font-medium text-slate-800">{{ tile.name }}</span>
            </label>
          </div>
        </div>
      </template>

      <template #footer>
        <div class="flex w-full justify-end gap-2 px-6 pb-6">
          <UButton
            variant="ghost"
            color="neutral"
            :label="t('inventory.dashboard_tiles.cancel')"
            :disabled="props.isSaving"
            @click="isOpen = false"
          />
          <UButton
            color="primary"
            :label="t('inventory.dashboard_tiles.save')"
            :loading="props.isSaving"
            :disabled="props.isSaving || !canSave"
            @click="saveSelection"
          />
        </div>
      </template>
    </UModal>
  </div>
</template>
