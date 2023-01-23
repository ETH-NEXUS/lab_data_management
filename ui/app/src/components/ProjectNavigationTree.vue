<script setup lang="ts">
import {ref, onMounted, computed, watchEffect} from 'vue'
import {handleError} from '../helpers/errorHandling'
import {Project, Experiment, Plate} from './models'
import {QTreeNode} from 'quasar'
import {useI18n} from 'vue-i18n'
import {useRouter} from 'vue-router'
import {useQuasar} from 'quasar'
import {useSettingsStore} from '../stores/settings'
import {storeToRefs} from 'pinia'
import {useProjectStore} from '../stores/project'

const router = useRouter()
const {t} = useI18n()
const $q = useQuasar()

const projectStore = useProjectStore()

const initialize = async () => {
  try {
    await projectStore.initialize()
    updateProjectNodes()
  } catch (err) {
    handleError(err)
  }
}

onMounted(async () => {
  initialize()
})

const {navigationTree, projectNavigationTree} = storeToRefs(useSettingsStore())

watchEffect(() => {
  if (projectNavigationTree.value.needsUpdate) {
    initialize()
    projectNavigationTree.value.needsUpdate = false
  }
})

const nodeHandler = (node: QTreeNode) => {
  if ('plate' in node) {
    router.push(`/plate/${node.plate.barcode}`)
  } else if ('project' in node) {
    router.push(`/project/${node.project.name}`)
  }
}

const projectNodes = ref<QTreeNode>({
  label: t('label.projects'),
  icon: 'biotech',
  header: 'projects',
  children: [],
})

const addProjectNode = (project: Project) => {
  const node: QTreeNode = {
    label: project.name,
    icon: '0_work',
    header: 'project',
    children: [],
    project: project,
  }
  projectNodes.value.children?.push(node)
  for (const experiment of project.experiments) {
    addExperimentNode(project, experiment)
    for (const plate of experiment.plates) {
      addPlateNode(experiment, plate)
    }
    sortPlateNodes(experiment)
  }
}

const addExperimentNode = (project: Project, experiment: Experiment) => {
  const projectNode = projectNodes.value.children?.find(c => c.project.id === project.id)
  if (projectNode) {
    projectNode.children?.push({
      label: `${experiment.name}`,
      icon: 'o_science',
      header: 'experiment',
      handler: nodeHandler,
      children: [],
      experiment: experiment,
    })
  } else {
    handleError(`TSNH: Project ${project.name} not found in tree.`)
  }
}

const addPlateNode = (experiment: Experiment, plate: Plate) => {
  const projectNode = projectNodes.value.children?.find(c => c.project.id === experiment.project)
  if (projectNode) {
    const experimentNode = projectNode.children?.find(c => c.experiment.id === experiment.id)
    if (experimentNode) {
      experimentNode.children?.push({
        label: `${plate.barcode} (${plate.dimension?.name || t('message.no_dimension')})`,
        icon: 'o_view_module',
        header: 'plate',
        handler: nodeHandler,
        plate: plate,
      })
    } else {
      handleError(`TSNH: Experiment ${experiment.name} not found in tree.`)
    }
  } else {
    handleError(`TSNH: Project ${experiment.project} not found in tree.`)
  }
}

const sortPlateNodes = (experiment: Experiment) => {
  const projectNode = projectNodes.value.children?.find(c => c.project.id === experiment.project)
  if (projectNode) {
    const experimentNode = projectNode.children?.find(c => c.experiment.id === experiment.id)
    if (experimentNode) {
      experimentNode.children = experimentNode.children?.sort((n1, n2) =>
        n1.plate.barcode.localeCompare(n2.plate.barcode)
      )
    }
  }
}

const updateProjectNodes = () => {
  projectNodes.value.children = []
  if (projectStore.projects) {
    for (const project of projectStore.projects) {
      addProjectNode(project)
    }
  }
}

const nodes = computed<Array<QTreeNode>>(() => {
  const nodes: Array<QTreeNode> = []
  nodes.push(projectNodes.value)
  return nodes
})

const newProject = async () => {
  $q.dialog({
    title: t('title.project_name'),
    message: t('message.project_name'),
    prompt: {
      model: '',
      type: 'text',
    },
    cancel: true,
    persistent: true,
  }).onOk(async projectName => {
    if (projectNodes.value.children) {
      try {
        const project = await projectStore.add(projectName)
        addProjectNode(project)
      } catch (err) {
        handleError(err)
      }
    }
  })
}

const newExperiment = async (project: Project) => {
  $q.dialog({
    title: t('title.experiment_name'),
    message: t('message.experiment_name'),
    prompt: {
      model: '',
      type: 'text',
    },
    cancel: true,
    persistent: true,
  }).onOk(async experimentName => {
    if (projectNodes.value.children) {
      try {
        const experiment = await projectStore.addExperiment(project, experimentName)
        addExperimentNode(project, experiment)
      } catch (err) {
        handleError(err, false)
      }
    }
  })
}

const newPlate = async (experiment: Experiment) => {
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
    if (projectNodes.value.children) {
      try {
        const plate = await projectStore.addPlate(experiment, barcode)
        addPlateNode(experiment, plate)
        sortPlateNodes(experiment)
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
    v-model:expanded="projectNavigationTree.expandedNodes"
    :filter="navigationTree.filter">
    <template v-slot:header-projects="prop">
      <q-icon :name="prop.node.icon || 'star'" size="28px" class="q-mr-sm" />
      <q-menu touch-position context-menu>
        <q-list dense style="min-width: 100px">
          <q-item clickable v-close-popup>
            <q-item-section @click="newProject">{{ t('action.new_project') }}</q-item-section>
          </q-item>
        </q-list>
      </q-menu>
      {{ prop.node.label }}
    </template>
    <template v-slot:header-project="prop">
      <q-icon :name="prop.node.icon || 'star'" size="24px" class="q-mr-sm" style="justify-content: end" />
      <q-menu touch-position context-menu>
        <q-list dense style="min-width: 100px">
          <q-item clickable v-close-popup>
            <q-item-section>{{ t('action.project_properties') }}</q-item-section>
          </q-item>
          <q-separator />
          <q-item clickable v-close-popup>
            <q-item-section @click="newExperiment(prop.node.project)">
              {{ t('action.new_experiment') }}
            </q-item-section>
          </q-item>
        </q-list>
      </q-menu>
      {{ prop.node.label }}
    </template>
    <template v-slot:header-experiment="prop">
      <q-icon :name="prop.node.icon || 'star'" size="24px" class="q-mr-sm" style="justify-content: end" />
      <q-menu touch-position context-menu>
        <q-list dense style="min-width: 100px">
          <q-item clickable v-close-popup>
            <q-item-section>{{ t('action.experiment_properties') }}</q-item-section>
          </q-item>
          <q-separator />
          <q-item clickable v-close-popup>
            <q-item-section @click="newPlate(prop.node.experiment)">
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
