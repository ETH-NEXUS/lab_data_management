<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import type { InventoryDashboardTilePreference } from '~/types/inventory'

type Props = {
  tiles: InventoryDashboardTilePreference[]
  isLoading: boolean
  isSaving: boolean
  hasError: boolean
  saveVersion: number
}

const props = defineProps<Props>()

const emit = defineEmits<{
  (e: 'save', tileKeys: string[]): void
}>()

const { t } = useI18n()
const isOpen = ref(false)
const selectedTileKeys = ref<string[]>([])

const setSelectedTileKeys = (tiles: InventoryDashboardTilePreference[]): void => {
  selectedTileKeys.value = tiles
    .filter((tile) => tile.is_visible)
    .sort((leftTile, rightTile) => leftTile.position - rightTile.position)
    .map((tile) => tile.key)
}

const selectedTiles = computed<InventoryDashboardTilePreference[]>(() => {
  return selectedTileKeys.value
    .map((tileKey) => props.tiles.find((tile) => tile.key === tileKey))
    .filter((tile): tile is InventoryDashboardTilePreference => tile !== undefined)
})

const availableTiles = computed<InventoryDashboardTilePreference[]>(() => {
  return props.tiles
    .filter((tile) => !isTileSelected(tile.key))
    .sort((leftTile, rightTile) => leftTile.position - rightTile.position)
})

watch(
  () => props.tiles,
  (tiles) => {
    if (!isOpen.value) {
      setSelectedTileKeys(tiles)
    }
  },
  { immediate: true },
)

watch(isOpen, (isModalOpen) => {
  if (isModalOpen) {
    setSelectedTileKeys(props.tiles)
  }
})

watch(
  () => props.saveVersion,
  () => {
    isOpen.value = false
  },
)

const isTileSelected = (tileKey: string): boolean => selectedTileKeys.value.includes(tileKey)

const toggleTile = (tileKey: string, isSelected: boolean): void => {
  if (isSelected && !isTileSelected(tileKey)) {
    selectedTileKeys.value.push(tileKey)
  }

  if (!isSelected && isTileSelected(tileKey)) {
    selectedTileKeys.value = selectedTileKeys.value.filter((selectedTileKey) => selectedTileKey !== tileKey)
  }
}

const moveTile = (tileKey: string, offset: -1 | 1): void => {
  const currentIndex = selectedTileKeys.value.indexOf(tileKey)
  const nextIndex = currentIndex + offset

  if (currentIndex < 0 || nextIndex < 0 || nextIndex >= selectedTileKeys.value.length) return

  const nextTileKeys = [...selectedTileKeys.value]
  const [tileToMove] = nextTileKeys.splice(currentIndex, 1)
  nextTileKeys.splice(nextIndex, 0, tileToMove)
  selectedTileKeys.value = nextTileKeys
}

const saveSelection = (): void => {
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
      :disabled="props.isLoading || props.hasError"
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
            {{ t('inventory.dashboard_tiles.selected_count', { count: selectedTileKeys.length }) }}
          </p>

          <div class="space-y-2">
            <p class="text-sm font-semibold text-slate-800">
              {{ t('inventory.dashboard_tiles.selected_title') }}
            </p>
            <div
              v-for="(tile, index) in selectedTiles"
              :key="tile.key"
              class="flex items-center gap-3 rounded-lg border border-slate-200 px-3 py-2"
            >
              <label class="flex min-w-0 flex-1 cursor-pointer items-center gap-3">
                <UCheckbox
                  :model-value="true"
                  :disabled="props.isSaving"
                  @update:model-value="(isSelected) => toggleTile(tile.key, Boolean(isSelected))"
                />
                <span class="truncate text-sm font-medium text-slate-800">{{ tile.name }}</span>
              </label>
              <UButton
                variant="ghost"
                color="neutral"
                icon="i-heroicons-arrow-up"
                :aria-label="t('inventory.dashboard_tiles.move_up')"
                :title="t('inventory.dashboard_tiles.move_up')"
                :disabled="props.isSaving || index === 0"
                @click="moveTile(tile.key, -1)"
              />
              <UButton
                variant="ghost"
                color="neutral"
                icon="i-heroicons-arrow-down"
                :aria-label="t('inventory.dashboard_tiles.move_down')"
                :title="t('inventory.dashboard_tiles.move_down')"
                :disabled="props.isSaving || index === selectedTiles.length - 1"
                @click="moveTile(tile.key, 1)"
              />
            </div>
            <p
              v-if="selectedTiles.length === 0"
              class="rounded-lg border border-dashed border-slate-200 px-3 py-2 text-sm text-slate-600"
            >
              {{ t('inventory.dashboard_tiles.no_selected_tiles') }}
            </p>
          </div>

          <div v-if="availableTiles.length > 0" class="space-y-2">
            <p class="text-sm font-semibold text-slate-800">
              {{ t('inventory.dashboard_tiles.available_title') }}
            </p>
            <label
              v-for="tile in availableTiles"
              :key="tile.key"
              class="flex cursor-pointer items-center gap-3 rounded-lg border border-slate-200 px-3 py-2 hover:bg-slate-50"
            >
              <UCheckbox
                :model-value="false"
                :disabled="props.isSaving"
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
            :disabled="props.isSaving"
            @click="saveSelection"
          />
        </div>
      </template>
    </UModal>
  </div>
</template>
