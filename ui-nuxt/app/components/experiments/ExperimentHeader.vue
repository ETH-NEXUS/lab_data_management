<script setup lang="ts">
import { computed } from 'vue'
import type { Experiment } from '~/types/experiments'

const props = defineProps<{
  experiment: Experiment
  projectName?: string | null
}>()

const emit = defineEmits<{
  (e: 'edit-name'): void
}>()

const { t } = useI18n()

const title = computed(() => {
  if (props.projectName) {
    return t('experiments.header.title_with_project', {
      projectName: props.projectName,
      experimentName: props.experiment.name,
    })
  }

  return t('experiments.header.title', { experimentName: props.experiment.name })
})
</script>

<template>
  <div class="flex items-center justify-center gap-2">
    <h1 class="text-primary text-center text-4xl font-semibold">{{ title }}</h1>

    <UButton
      size="xs"
      variant="ghost"
      icon="i-heroicons-pencil-square"
      :aria-label="t('experiments.header.edit_name_aria')"
      @click="emit('edit-name')"
    />
  </div>
</template>
