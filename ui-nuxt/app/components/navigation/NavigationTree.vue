<script setup lang="ts">
import { ref } from 'vue'
import LibraryNavigationTree from '~/components/navigation/LibraryNavigationTree.vue'
import ManagementNavigationTree from '~/components/navigation/ManagementNavigationTree.vue'
import ProjectNavigationTree from '~/components/navigation/ProjectNavigationTree.vue'
import TemplateNavigationTree from '~/components/navigation/TemplateNavigationTree.vue'

const filter = ref('')
const includeDeprecatedFunctionality = ref(true)

const resetFilter = (): void => {
  filter.value = ''
}
</script>

<template>
  <div class="p-4 space-y-3">
    <UInput
      v-model="filter"
      placeholder="Filter navigation"
      icon="i-heroicons-magnifying-glass"
    >
      <template #trailing>
        <UButton
          v-if="filter !== ''"
          size="xs"
          variant="ghost"
          icon="i-heroicons-x-mark"
          aria-label="Clear filter"
          @click="resetFilter"
        />
      </template>
    </UInput>

    <ProjectNavigationTree :filter="filter" />
    <LibraryNavigationTree :filter="filter" />
    <TemplateNavigationTree v-if="includeDeprecatedFunctionality" :filter="filter" />
    <ManagementNavigationTree :filter="filter" />
  </div>
</template>
