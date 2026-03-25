<script setup lang="ts">
import {ref} from 'vue'
import {useI18n} from 'vue-i18n'
import {useSettingsStore} from 'stores/settings'
import {storeToRefs} from 'pinia'
import LibraryNavigationTree from './LibraryNavigationTree.vue'
import ProjectNavigationTree from './ProjectNavigationTree.vue'
import TemplateNavigationTree from './TemplateNavigationTree.vue'
import ManagementNavigationTree from 'components/management/ManagementNavigationTree.vue'

const {t} = useI18n()

const {navigationTree, include_depricated_functionality} = storeToRefs(useSettingsStore())
const filterRef = ref<HTMLInputElement>()

const resetFilter = () => {
  navigationTree.value.filter = ''
  filterRef.value?.focus()
}
</script>

<template>
  <div class="q-pa-md q-gutter-sm">
    <q-input ref="filterRef" v-model="navigationTree.filter" :label="t('label.filter')">
      <template v-slot:append>
        <q-icon
          v-if="navigationTree.filter !== ''"
          name="clear"
          class="cursor-pointer"
          @click="resetFilter" />
      </template>
    </q-input>

    <project-navigation-tree />
    <library-navigation-tree />
    <template-navigation-tree v-if="include_depricated_functionality" />
    <management-navigation-tree />
  </div>
</template>

<style lang="sass">
.q-tree__node-header-content
  cursor: pointer
.link:hover
  color: $primary
</style>
