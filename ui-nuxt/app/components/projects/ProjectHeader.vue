<script setup lang="ts">
import { computed } from 'vue'
import type { Project } from '~/types/projects'

const props = defineProps<{
  project: Project
}>()

const emit = defineEmits<{
  (e: 'edit-name'): void
}>()

const { t } = useI18n()

const title = computed(() => t('projects.header.title', { name: props.project.name }))
</script>

<template>
  <div class="flex items-center justify-center gap-2">
    <h1 class="text-primary text-center text-4xl font-semibold">{{ title }}</h1>

    <UButton
      v-if="!props.project.harvest_id"
      size="xs"
      variant="ghost"
      icon="i-heroicons-pencil-square"
      :aria-label="t('projects.header.edit_name_aria')"
      @click="emit('edit-name')"
    />
  </div>
</template>
