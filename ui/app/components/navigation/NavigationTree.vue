<script setup lang="ts">
import { ref } from 'vue'
import InventoryNavigationTree from '~/components/navigation/InventoryNavigationTree.vue'
import LibraryNavigationTree from '~/components/navigation/LibraryNavigationTree.vue'
import ManagementNavigationTree from '~/components/navigation/ManagementNavigationTree.vue'
import ProjectNavigationTree from '~/components/navigation/ProjectNavigationTree.vue'

const filter = ref('')

const resetFilter = (): void => {
  filter.value = ''
}
</script>

<template>
  <div class="space-y-3 p-4">
    <UInput v-model="filter" placeholder="Filter navigation" icon="i-heroicons-magnifying-glass">
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
    <InventoryNavigationTree :filter="filter" />
    <LibraryNavigationTree :filter="filter" />
    <ManagementNavigationTree :filter="filter" />
  </div>
</template>
