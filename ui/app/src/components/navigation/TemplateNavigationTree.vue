<script setup lang="ts">
import {ref, onMounted, computed, watchEffect} from 'vue'
import {handleError} from '../../helpers/errorHandling'
import {TemplateCategory, Template} from '../models'
import {QTreeNode} from 'quasar'
import {useI18n} from 'vue-i18n'
import {useRouter} from 'vue-router'
import {useQuasar} from 'quasar'
import {useSettingsStore} from 'stores/settings'
import {storeToRefs} from 'pinia'
import {useTemplateStore} from 'stores/template'

const router = useRouter()
const {t} = useI18n()
const $q = useQuasar()

const templateStore = useTemplateStore()

const initialize = async () => {
  try {
    await templateStore.initialize()
    updateCategoryNodes()
  } catch (err) {
    handleError(err)
  }
}

onMounted(async () => {
  initialize()
})

const {navigationTree, templateNavigationTree} = storeToRefs(useSettingsStore())

watchEffect(() => {
  if (templateNavigationTree.value.needsUpdate) {
    initialize()
    templateNavigationTree.value.needsUpdate = false
  }
})

const nodeHandler = async (node: QTreeNode) => {
  if ('template' in node) {
    if (!node.template.plate) {
      // for convenience reasons we add a plate to the template if it's missing
      const plate = await templateStore.addPlate(node.template, node.category.name)
      node.template.plate = plate
    }
    router.push(`/plate/${node.template.plate.barcode}`)
  }
}

const categoryNodes = ref<QTreeNode>({
  label: t('label.templates'),
  icon: 'category',
  header: 'categories',
  children: [],
})

const addCategoryNode = (category: TemplateCategory) => {
  const node: QTreeNode = {
    label: category.name,
    icon: 'o_hexagon',
    header: 'category',
    children: [],
    category: category,
  }
  categoryNodes.value.children?.push(node)
  for (const template of category.templates) {
    addTemplateNode(category, template)
    sortTemplateNodes(category)
  }
}

const addTemplateNode = (category: TemplateCategory, template: Template) => {
  const categoryNode = categoryNodes.value.children?.find(c => c.category.id === category.id)
  if (categoryNode) {
    categoryNode.children?.push({
      label: `${template.name} (${template.plate.dimension || t('message.no_dimension')})`,
      icon: 'o_view_module',
      header: 'template',
      handler: nodeHandler,
      children: [],
      category: category,
      template: template,
    })
  } else {
    handleError(`TSNH: Category ${category.name} not found in tree.`)
  }
}

const sortTemplateNodes = (category: TemplateCategory) => {
  const categoryNode = categoryNodes.value.children?.find(c => c.category.id === category.id)
  if (categoryNode) {
    categoryNode.children = categoryNode.children?.sort((n1, n2) =>
      n1.template.name.localeCompare(n2.template.name)
    )
  }
}

const updateCategoryNodes = () => {
  categoryNodes.value.children = []
  if (templateStore.templateCategories) {
    for (const category of templateStore.templateCategories) {
      addCategoryNode(category)
    }
  }
}

const nodes = computed<Array<QTreeNode>>(() => {
  const nodes: Array<QTreeNode> = []
  nodes.push(categoryNodes.value)
  return nodes
})

const newCategory = async () => {
  $q.dialog({
    title: t('title.project_name'),
    message: t('message.project_name'),
    prompt: {
      model: '',
      type: 'text',
    },
    cancel: true,
    persistent: true,
  }).onOk(async categoryName => {
    if (categoryNodes.value.children) {
      try {
        const project = await templateStore.addCategory(categoryName)
        addCategoryNode(project)
      } catch (err) {
        handleError(err)
      }
    }
  })
}

const newTemplate = async (category: TemplateCategory) => {
  $q.dialog({
    title: t('title.template_name'),
    message: t('message.template_name'),
    prompt: {
      model: '',
      type: 'text',
    },
    cancel: true,
    persistent: true,
  }).onOk(async templateName => {
    if (categoryNodes.value.children) {
      try {
        const template = await templateStore.addTemplate(category, templateName)
        addTemplateNode(category, template)
      } catch (err) {
        handleError(err, false)
      }
    }
  })
}
</script>

<template>
  <q-tree
    :nodes="nodes"
    dense
    node-key="label"
    v-model:expanded="templateNavigationTree.expandedNodes"
    :filter="navigationTree.filter">
    <template v-slot:header-categories="prop">
      <q-icon :name="prop.node.icon || 'star'" size="28px" class="q-mr-sm" />
      <q-menu touch-position context-menu>
        <q-list dense style="min-width: 100px">
          <q-item clickable v-close-popup>
            <q-item-section @click="newCategory">{{ t('action.new_category') }}</q-item-section>
          </q-item>
        </q-list>
      </q-menu>
      {{ prop.node.label }}
    </template>
    <template v-slot:header-category="prop">
      <q-icon :name="prop.node.icon || 'star'" size="24px" class="q-mr-sm" style="justify-content: end" />
      <q-menu touch-position context-menu>
        <q-list dense style="min-width: 100px">
          <q-item clickable v-close-popup>
            <q-item-section>{{ t('action.category_properties') }}</q-item-section>
          </q-item>
          <q-separator />
          <q-item clickable v-close-popup>
            <q-item-section @click="newTemplate(prop.node.category)">
              {{ t('action.new_template') }}
            </q-item-section>
          </q-item>
        </q-list>
      </q-menu>
      {{ prop.node.label }}
    </template>
  </q-tree>
</template>

<style lang="sass">
.q-tree__node-header-content
  cursor: pointer
</style>
