<script setup lang="ts">
import BaseButton from '~/components/common/BaseButton.vue'
import type { Threshold } from '~/types/messages'

const props = defineProps<{
  threshold: Threshold
  isRecalculating?: boolean
}>()

const emit = defineEmits<{
  (e: 'edit-threshold' | 'recalculate-status'): void
}>()

const { t } = useI18n()

const onEditThreshold = () => emit('edit-threshold')
const onRecalculateStatus = () => emit('recalculate-status')
</script>

<template>
  <UCard
    class="mx-auto w-[80%]"
    :ui="{
      root: 'border border-white/40 bg-white/30 backdrop-blur-md shadow-sm divide-y divide-white/20',
    }"
  >
    <template #header>
      <p class="text-primary font-semibold">{{ t('messages_page.sections.thresholds.title') }}</p>
    </template>

    <div class="space-y-4">
      <div class="space-y-1 text-sm text-slate-700">
        <p>
          <span class="font-semibold">{{ t('messages_page.sections.thresholds.dmso') }}:</span>
          {{ props.threshold.dmso }}%
        </p>
        <p>
          <span class="font-semibold">{{ t('messages_page.sections.thresholds.amount') }}:</span>
          {{ props.threshold.amount }}µL
        </p>
      </div>

      <div class="flex flex-wrap gap-2">
        <BaseButton
          :label="t('messages_page.sections.thresholds.edit_button')"
          :on-click="onEditThreshold"
          variant="secondary"
          size="sm"
          width="auto"
          class-name="text-teal-800 hover:text-teal-900 hover:bg-teal-50 hover:border-teal-200"
        />

        <BaseButton
          :label="t('messages_page.sections.thresholds.recalculate_button')"
          :on-click="onRecalculateStatus"
          variant="primary"
          size="sm"
          width="auto"
          :loading="props.isRecalculating"
        />
      </div>
    </div>
  </UCard>
</template>
