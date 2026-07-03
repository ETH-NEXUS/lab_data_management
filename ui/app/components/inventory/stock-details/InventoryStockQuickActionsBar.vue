<script setup lang="ts">
type Props = {
  isFavorite: boolean
  isEditingStock: boolean
  isMovingStock: boolean
  isRecordingUsage: boolean
  isTogglingFavorite: boolean
  isArchivingStock: boolean
}

const props = defineProps<Props>()

const emit = defineEmits<{
  (e: 'adjust-stock' | 'move-item' | 'record-usage' | 'toggle-favorite' | 'archive-item'): void
}>()

const { t } = useI18n()
</script>

<template>
  <div class="flex flex-col gap-3 sm:flex-row sm:flex-wrap sm:items-center sm:justify-between">
    <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
      {{ t('inventory.page.section_labels.quick_actions') }}
    </p>
    <div class="flex flex-wrap items-center gap-2">
      <UButton
        color="primary"
        size="xs"
        :label="t('inventory.stock_drawer.actions.adjust_stock')"
        :disabled="props.isEditingStock"
        @click="emit('adjust-stock')"
      />
      <UButton
        color="neutral"
        variant="soft"
        size="xs"
        :label="t('inventory.stock_drawer.actions.move_item')"
        :disabled="props.isMovingStock"
        @click="emit('move-item')"
      />
      <UButton
        color="neutral"
        variant="soft"
        size="xs"
        label="Record usage"
        :disabled="props.isRecordingUsage"
        @click="emit('record-usage')"
      />
      <UButton
        color="warning"
        variant="soft"
        size="xs"
        :label="props.isFavorite ? 'Unfavorite item' : 'Favorite item'"
        :disabled="props.isTogglingFavorite"
        :loading="props.isTogglingFavorite"
        @click="emit('toggle-favorite')"
      />
      <UButton
        color="neutral"
        variant="soft"
        size="xs"
        label="Archive item"
        :disabled="props.isArchivingStock"
        :loading="props.isArchivingStock"
        @click="emit('archive-item')"
      />
    </div>
  </div>
</template>
