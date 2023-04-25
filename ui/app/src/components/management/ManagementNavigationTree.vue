<script setup lang="ts">
import {useManagementStore} from 'stores/management'
import {onMounted, computed, ref} from 'vue'
import {storeToRefs} from 'pinia'
import {useRouter} from 'vue-router'
import {FileSystemItem} from 'components/models'
import {QTreeNode, useQuasar} from 'quasar'
import {useI18n} from 'vue-i18n'

const {t} = useI18n()
const $q = useQuasar()

const {dataDirectory, selectedPath, selectedPaths} = storeToRefs(useManagementStore())
const managementStore = useManagementStore()
const router = useRouter()
const deleteDialog = ref<boolean>(false)

const nodes = computed<QTreeNode[]>(() => {
  const res = convertToQTreeNodes([dataDirectory.value])

  return [{handler: nodeHandler, key: 'Management', label: 'Management', icon: 'troubleshoot', children: res}]
})

onMounted(async () => {
  await managementStore.initialize()
})

const treeLabel = (node: FileSystemItem) => {
  return node.name
}

const treeChildren = (node: FileSystemItem) => {
  if (node.type === 'directory') {
    return node.children
  }
}

const convertToQTreeNodes = (fileSystemItems: FileSystemItem[]): QTreeNode[] => {
  return fileSystemItems.map(item => {
    const qTreeNode: QTreeNode = {
      handler: nodeHandler,
      header: 'paths',
      icon: item.type === 'directory' ? 'folder' : 'description',
      key: item.name,
      label: item.name,
      type: item.type,
      name: item.name,
      path: item.path,
      children: item.type === 'directory' ? convertToQTreeNodes(item.children) : [],
    }
    return qTreeNode
  })
}

const nodeHandler = async (node: QTreeNode) => {
  if (node.label === 'Management') {
    await router.push('/management')
  } else if (node.type === 'file' || node.type === 'directory') {
    selectedPaths.value.push(node.path)
  }
}

const handleDeleteDialog = (path: string) => {
  selectedPath.value = path
  deleteDialog.value = true
}

const deleteFile = async () => {
  await managementStore.deleteFile(selectedPath.value)
  deleteDialog.value = false
}
</script>

<template>
  <q-tree
    v-model:selected="selectedPath"
    selected-color="primary"
    animated
    dense
    v-if="nodes"
    :nodes="nodes"
    node-key="name"
    :label-fn="treeLabel"
    :children-fn="treeChildren">
    <template v-slot:header-paths="prop">
      <q-icon :name="prop.node.icon" size="20px" class="q-mr-sm" />
      <q-menu touch-position context-menu>
        <q-list dense style="min-width: 100px">
          <q-item clickable v-close-popup>
            <q-item-section @click="handleDeleteDialog(prop.node.path)">
              {{ t('action.delete_item') }}
            </q-item-section>
          </q-item>
        </q-list>
      </q-menu>
      <span class="fileSystemItem">
        {{ prop.node.label }}
      </span>
    </template>
  </q-tree>
  <q-dialog v-model="deleteDialog">
    <q-card>
      <q-card-section>
        <div>Are you sure you want to delete the item {{ selectedPath }}?</div>
      </q-card-section>
      <q-card-actions align="right">
        <q-btn flat label="OK" color="primary" v-close-popup @click="deleteFile"></q-btn>
      </q-card-actions>
    </q-card>
  </q-dialog>
</template>

<style lang="sass">
.q-tree__node-header-content
  cursor: pointer


.fileSystemItem:hover
    color: $primary
</style>
