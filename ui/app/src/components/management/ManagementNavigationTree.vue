<script setup lang="ts">
import {useManagementStore} from 'stores/management'
import {onMounted, computed} from 'vue'
import {storeToRefs} from 'pinia'
import {useRouter} from 'vue-router'
import {FileSystemItem} from 'components/models'
import {QTreeNode} from 'quasar'

const {dataDirectory, selectedPath, selectedPaths} = storeToRefs(useManagementStore())
const managementStore = useManagementStore()
const router = useRouter()

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
    <template v-slot:header-categories="prop">
      <q-icon name="troubleshoot" size="28px" class="q-mr-sm" />
      <span class="fileSystemItem">
        {{ prop.node.label }}
      </span>
    </template>
  </q-tree>
</template>

<style lang="sass">
.q-tree__node-header-content
  cursor: pointer
.fileSystemItem
  color: red

.fileSystemItem:hover
    color: $primary
</style>
