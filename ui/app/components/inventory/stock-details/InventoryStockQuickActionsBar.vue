<script setup lang="ts">
type Props = {
  isEditingStock: boolean
  isMovingStock: boolean
  isRecordingUsage: boolean
}

const props = defineProps<Props>()

const emit = defineEmits<{
  (e: 'adjust-stock' | 'move-item' | 'record-usage'): void
}>()

const { t } = useI18n()
</script>

<template>
  <div class="flex items-center justify-between gap-3">
    <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
      {{ t('inventory.page.section_labels.quick_actions') }}
    </p>
    <div class="flex items-center gap-2">
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
    </div>
  </div>
</template>
