<script setup lang="ts">
import { storeToRefs } from 'pinia'
import { useManagementStore } from '~/stores/management'

const { t } = useI18n()
const managementStore = useManagementStore()
const { selectedPaths } = storeToRefs(managementStore)

/**
 * Removes one selected path from the list.
 *
 * Accepted data example:
 * - `'/data/imports/file.csv'`
 */
const deletePath = (path: string): void => {
  managementStore.removeSelectedPath(path)
}

/**
 * Adds drag payload with one selected path.
 *
 * Accepted data example:
 * - `path = '/data/imports/file.csv'`
 */
const handleDragStart = (path: string, event: DragEvent): void => {
  const normalizedPath = path.trim()
  if (normalizedPath === '' || !event.dataTransfer) {
    return
  }

  // Keep drag payload compatible across desktop browsers and drop targets.
  event.dataTransfer.effectAllowed = 'copy'
  event.dataTransfer.setData('application/x-management-path', normalizedPath)
  event.dataTransfer.setData('text/plain', normalizedPath)
  event.dataTransfer.setData('text', normalizedPath)
}
</script>

<template>
  <UCard
    class="mx-auto w-full max-w-5xl"
    :ui="{
      root: 'core-card divide-y divide-slate-200/70',
    }"
  >
    <template #header>
      <p class="font-semibold text-blue-700">{{ t('management.selected_directories') }}</p>
    </template>

    <p v-if="selectedPaths.length === 0" class="text-sm text-slate-600">
      {{ t('management.no_dirs_selected') }}
    </p>

    <ul v-else class="space-y-2">
      <li
        v-for="(path, index) in selectedPaths"
        :key="`${path}-${index}`"
        class="flex items-start justify-between gap-3 rounded-md border border-slate-200 bg-white px-3 py-2"
      >
        <span draggable="true" class="min-w-0 text-sm text-slate-700" @dragstart="handleDragStart(path, $event)">
          <b class="mr-2">{{ `${index + 1}.` }}</b>
          <i class="break-all">{{ path }}</i>
        </span>

        <button
          type="button"
          class="inline-flex h-7 w-7 shrink-0 cursor-pointer items-center justify-center rounded-md border border-blue-100 bg-blue-50 text-blue-700 transition hover:bg-blue-100"
          :aria-label="t('management.remove_selected_path_aria')"
          @click="deletePath(path)"
        >
          <UIcon name="i-heroicons-trash" class="h-4 w-4" />
        </button>
      </li>
    </ul>
  </UCard>
</template>
