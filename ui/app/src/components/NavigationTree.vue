<script setup lang="ts">
import {ref, onMounted, computed} from 'vue'
import {api} from '../boot/axios'
import {handleError} from '../helpers/errorHandling'
import {CompoundLibrary, Project, Experiment, Plate} from './models'
import {QTreeNode} from 'quasar'
import {useI18n} from 'vue-i18n'
import {useRouter} from 'vue-router'
import {useQuasar} from 'quasar'
import {useSettingsStore} from '../stores/settings'
import {storeToRefs} from 'pinia'

const router = useRouter()
const {t} = useI18n()
const $q = useQuasar()

const compoundLibraries = ref<Array<CompoundLibrary> | undefined>()
const projects = ref<Array<Project> | undefined>()

onMounted(async () => {
  try {
    const resp_cl = await api.get('/api/compoundlibraries/')
    compoundLibraries.value = resp_cl.data.results
    updateCompoundLibraryNodes()

    const resp_p = await api.get('/api/projects/')
    projects.value = resp_p.data.results
    updateProjectNodes()
  } catch (err) {
    handleError(err)
  }
})

const {navigationTree} = storeToRefs(useSettingsStore())
const filter = ref('')
const filterRef = ref<HTMLInputElement>()

const resetFilter = () => {
  filter.value = ''
  filterRef.value?.focus()
}

const nodeHandler = (node: QTreeNode) => {
  if ('plate' in node) {
    router.push(`/plate/${node.plate.barcode}`)
  } else if ('project' in node) {
    router.push(`/project/${node.project.name}`)
  }
}

const compoundLibraryNodes = ref<QTreeNode>({
  label: t('label.compound_libraries'),
  icon: 'science',
  header: 'compoundlibs',
  children: [],
})
const projectNodes = ref<QTreeNode>({
  label: t('label.projects'),
  icon: 'biotech',
  header: 'projects',
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
      label: `${plate.barcode} (${plate.dimension})`,
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
  if (compoundLibraries.value) {
    for (const library of compoundLibraries.value) {
      addCompoundLibraryNode(library)
    }
  }
}

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
        label: `${plate.barcode}`,
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

const updateProjectNodes = () => {
  projectNodes.value.children = []
  if (projects.value) {
    for (const project of projects.value) {
      addProjectNode(project)
    }
  }
}

const nodes = computed<Array<QTreeNode>>(() => {
  const nodes: Array<QTreeNode> = []
  nodes.push(compoundLibraryNodes.value)
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
  })
    .onOk(async projectName => {
      if (projectNodes.value.children) {
        try {
          const resp = await api.post('/api/projects/', {
            name: projectName,
          })
          const project = resp.data
          projects.value?.push(project)
          addProjectNode(project)
        } catch (err) {
          handleError(err)
        }
      }
    })
    .onCancel(() => {
      // console.log('>>>> Cancel')
    })
    .onDismiss(() => {
      // console.log('I am triggered on both OK and Cancel')
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
  })
    .onOk(async experimentName => {
      if (projectNodes.value.children) {
        try {
          const resp = await api.post('/api/experiments/', {
            name: experimentName,
            project: project.id,
          })
          const experiment = resp.data
          project.experiments.push(experiment)
          addExperimentNode(project, experiment)
        } catch (err) {
          handleError(err, false)
        }
      }
    })
    .onCancel(() => {
      // console.log('>>>> Cancel')
    })
    .onDismiss(() => {
      // console.log('I am triggered on both OK and Cancel')
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
  })
    .onOk(async barcode => {
      if (projectNodes.value.children) {
        try {
          const resp = await api.post('/api/plates/', {
            barcode: barcode,
            experiment: experiment.id,
          })
          const plate = resp.data
          experiment.plates.push(plate)
          addPlateNode(experiment, plate)
        } catch (err) {
          handleError(err, false)
        }
      }
    })
    .onCancel(() => {
      // console.log('>>>> Cancel')
    })
    .onDismiss(() => {
      // console.log('I am triggered on both OK and Cancel')
    })
}
</script>

<template>
  <div class="q-pa-md q-gutter-sm">
    <q-input ref="filterRef" v-model="filter" :label="t('label.filter')">
      <template v-slot:append>
        <q-icon v-if="filter !== ''" name="clear" class="cursor-pointer" @click="resetFilter" />
      </template>
    </q-input>
    <q-tree
      :nodes="nodes"
      dense
      node-key="label"
      v-model:expanded="navigationTree.expandedNodes"
      :filter="filter">
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
  </div>
</template>

<style lang="sass">
.q-tree__node-header-content
  cursor: pointer
</style>
