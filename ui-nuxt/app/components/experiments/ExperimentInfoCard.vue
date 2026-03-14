<script setup lang="ts">
import { computed } from 'vue'
import type { Experiment } from '~/types/experiments'

const props = defineProps<{
  experiment: Experiment
  createdAtLabel: string
}>()

const emit = defineEmits<{
  (e: 'edit-description'): void
}>()

const { t } = useI18n()

/**
 * Computes number of plates shown in the experiment info section.
 *
 * Returned data examples:
 * - `0`
 * - `12`
 */
const platesCount = computed(() => props.experiment.plates?.length ?? 0)
</script>

<template>
  <UCard class="mx-auto w-[80%]">
    <div class="space-y-4">
      <div class="max-w-[600px] break-words">
        <p class="text-sm font-semibold text-slate-700">{{ t('experiments.info.description') }}</p>
        <p class="text-sm text-slate-600">
          {{ props.experiment.description || t('experiments.info.no_description') }}
        </p>

        <UButton
          size="xs"
          variant="ghost"
          icon="i-heroicons-pencil-square"
          :aria-label="t('experiments.info.edit_description_aria')"
          @click="emit('edit-description')"
        />
      </div>

      <div>
        <p class="text-sm font-semibold text-slate-700">{{ t('experiments.info.number_of_plates') }}</p>
        <p class="text-sm text-slate-600">
          {{ platesCount }}
        </p>
      </div>

      <div>
        <p class="text-sm font-semibold text-slate-700">{{ t('experiments.info.created_at') }}</p>
        <p class="text-sm text-slate-600">
          {{ props.createdAtLabel }}
        </p>
      </div>
    </div>
  </UCard>
</template>
