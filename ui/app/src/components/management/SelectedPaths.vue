<script setup lang="ts">
import {storeToRefs} from 'pinia'
import {useManagementStore} from 'stores/management'
import {useI18n} from 'vue-i18n'

const {selectedPaths} = storeToRefs(useManagementStore())
const {t} = useI18n()

const deletePath = (path: string) => {
  const index = selectedPaths.value.indexOf(path)
  if (index > -1) {
    selectedPaths.value.splice(index, 1)
  }
}

const handleDragStart = (path: string, event: DragEvent) => {
  event.dataTransfer?.setData('text/plain', path)
}
</script>

<template>
  <div style="max-width: 750px" class="q-mb-lg">
    <div class="text-body1 q-mb-md">{{ t('management.selected_directories') }}:</div>
    <q-list dense padding class="rounded-borders" v-if="selectedPaths.length > 0">
      <q-item v-for="(path, index) in selectedPaths" :key="path + index">
        <q-item-section>
          <span :draggable="true" @dragstart="handleDragStart(path, $event)">
            <b>{{ `${index + 1}.&nbsp;&nbsp;` }}</b>
            <i>{{ `${path}` }}</i>
          </span>
        </q-item-section>
        <q-item-section avatar class="cursor-pointer" @click="deletePath(path)">
          <q-icon color="primary" name="delete"></q-icon>
        </q-item-section>
      </q-item>
    </q-list>
    <div v-else class="text-caption">{{ t('management.no_dirs_selected') }}</div>
  </div>
</template>

<style scoped></style>
