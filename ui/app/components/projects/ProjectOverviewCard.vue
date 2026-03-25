<script setup lang="ts">
import BaseButton from '~/components/common/BaseButton.vue'
import type { Project } from '~/types/projects'

const props = defineProps<{
  project: Project
  experimentsCount: number
  createdAtLabel: string
  harvestLoading?: boolean
}>()

const emit = defineEmits<{
  (e: 'edit-description' | 'update-harvest'): void
}>()

const { t } = useI18n()

const onUpdateHarvest = () => emit('update-harvest')
</script>

<template>
  <UCard
    class="mx-auto w-[80%]"
    :ui="{
      root: 'core-card divide-y divide-slate-200/70',
    }"
  >
    <div class="space-y-4">
      <div class="max-w-[600px] break-words">
        <p class="text-sm font-semibold text-slate-700">{{ t('projects.overview.description') }}</p>
        <p class="text-sm text-slate-600">
          {{ props.project.description || t('projects.overview.no_description') }}
        </p>

        <UButton
          size="xs"
          variant="ghost"
          icon="i-heroicons-pencil-square"
          :aria-label="t('projects.overview.edit_description_aria')"
          @click="emit('edit-description')"
        />
      </div>

      <div v-if="props.project.harvest_notes" class="max-w-[600px] break-words">
        <p class="text-sm font-semibold text-slate-700">{{ t('projects.overview.harvest_notes') }}</p>
        <p class="text-sm text-slate-600">
          {{ props.project.harvest_notes }}
        </p>
      </div>

      <div>
        <p class="text-sm font-semibold text-slate-700">{{ t('projects.overview.number_of_experiments') }}</p>
        <p class="text-sm text-slate-600">
          {{ props.experimentsCount }}
        </p>
      </div>

      <div>
        <p class="text-sm font-semibold text-slate-700">{{ t('projects.overview.created_at') }}</p>
        <p class="text-sm text-slate-600">
          {{ props.createdAtLabel }}
        </p>
      </div>

      <BaseButton
        v-if="props.project.harvest_id"
        :label="t('projects.overview.update_harvest')"
        :on-click="onUpdateHarvest"
        variant="primary"
        size="sm"
        width="auto"
        :loading="props.harvestLoading"
      />
    </div>
  </UCard>
</template>
