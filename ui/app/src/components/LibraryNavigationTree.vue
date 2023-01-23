<script setup lang="ts">
import {ref, onMounted, computed, watchEffect} from 'vue'
import {handleError} from '../helpers/errorHandling'
import {CompoundLibrary, Plate} from './models'
import {QTreeNode} from 'quasar'
import {useI18n} from 'vue-i18n'
import {useRouter} from 'vue-router'
import {useQuasar} from 'quasar'
import {useSettingsStore} from '../stores/settings'
import {storeToRefs} from 'pinia'
import {useCompoundLibraryStore} from '../stores/compoundLibrary'

const router = useRouter()
const {t} = useI18n()
const $q = useQuasar()

const compoundLibraryStore = useCompoundLibraryStore()

const initialize = async () => {
  try {
    await compoundLibraryStore.initialize()
    updateCompoundLibraryNodes()
  } catch (err) {
    handleError(err)
  }
}

onMounted(async () => {
  initialize()
})

const {navigationTree, libraryNavigationTree} = storeToRefs(useSettingsStore())

watchEffect(() => {
  if (libraryNavigationTree.value.needsUpdate) {
    initialize()
    libraryNavigationTree.value.needsUpdate = false
  }
})

const nodeHandler = (node: QTreeNode) => {
  if ('plate' in node) {
    router.push(`/plate/${node.plate.barcode}`)
  }
}

const compoundLibraryNodes = ref<QTreeNode>({
  label: t('label.compound_libraries'),
  icon: 'science',
  header: 'compoundlibs',
  children: [],
})

const addCompoundLibraryNode = (library: CompoundLibrary) => {
  const node: QTreeNode = {
    label: library.name,
    icon: 'science',
    header: 'compoundlib',
    children: [],
    library: library,
  }
  for (const plate of library.plates) {
    node.children?.push({
      label: `${plate.barcode} (${plate.dimension || t('message.no_dimension')}})`,
      icon: 'o_view_module',
      header: 'plate',
      handler: nodeHandler,
      plate: plate,
    })
  }
  compoundLibraryNodes.value.children?.push(node)
}

const updateCompoundLibraryNodes = () => {
  compoundLibraryNodes.value.children = []
  if (compoundLibraryStore.libraries) {
    for (const library of compoundLibraryStore.libraries) {
      addCompoundLibraryNode(library)
    }
  }
}

const addCompoundLibraryPlateNode = (library: CompoundLibrary, plate: Plate) => {
  const compoundLibraryNode = compoundLibraryNodes.value.children?.find(n => n.library.id === library.id)
  if (compoundLibraryNode) {
    compoundLibraryNode.children?.push({
      label: `${plate.barcode} (${plate.dimension})`,
      icon: 'o_view_module',
      header: 'plate',
      handler: nodeHandler,
      plate: plate,
    })
  } else {
    handleError(`TSNH: Library ${library.name} not found in tree.`)
  }
}

const sortCompoundLibraryPlateNodes = (library: CompoundLibrary) => {
  const compoundLibraryNode = compoundLibraryNodes.value.children?.find(n => n.library.id === library.id)
  if (compoundLibraryNode) {
    compoundLibraryNode.children = compoundLibraryNode.children?.sort((n1, n2) =>
      n1.plate.barcode.localeCompare(n2.plate.barcode)
    )
  }
}

const nodes = computed<Array<QTreeNode>>(() => {
  const nodes: Array<QTreeNode> = []
  nodes.push(compoundLibraryNodes.value)
  return nodes
})

const newCompoundLibPlate = async (library: CompoundLibrary) => {
  $q.dialog({
    title: t('title.plate_barcode'),
    message: t('message.plate_barcode'),
    prompt: {
      model: '',
      type: 'text',
    },
    cancel: true,
    persistent: true,
  }).onOk(async barcode => {
    try {
      const plate = await compoundLibraryStore.addPlate(library, barcode)
      addCompoundLibraryPlateNode(library, plate)
      sortCompoundLibraryPlateNodes(library)
    } catch (err) {
      handleError(err, false)
    }
  })
}
</script>

<template>
  <q-tree
    :nodes="nodes"
    dense
    node-key="label"
    v-model:expanded="libraryNavigationTree.expandedNodes"
    :filter="navigationTree.filter">
    <template v-slot:header-compoundlib="prop">
      <q-icon :name="prop.node.icon || 'star'" size="24px" class="q-mr-sm" style="justify-content: end" />
      <q-menu touch-position context-menu>
        <q-list dense style="min-width: 100px">
          <q-item clickable v-close-popup>
            <q-item-section>{{ t('action.compoundlib_properties') }}</q-item-section>
          </q-item>
          <q-separator />
          <q-item clickable v-close-popup>
            <q-item-section @click="newCompoundLibPlate(prop.node.library)">
              {{ t('action.new_plate') }}
            </q-item-section>
          </q-item>
        </q-list>
      </q-menu>
      {{ prop.node.label }}
    </template>
    <template v-slot:header-plate="prop">
      <q-icon :name="prop.node.icon || 'star'" size="24px" class="q-mr-sm" style="justify-content: end" />
      {{ prop.node.label }}
    </template>
  </q-tree>
</template>

<style lang="sass">
.q-tree__node-header-content
  cursor: pointer
</style>
