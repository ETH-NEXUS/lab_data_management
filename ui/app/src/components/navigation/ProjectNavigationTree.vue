<script setup lang="ts">
import {ref, onMounted, computed, watchEffect} from 'vue'
import {handleError} from '../../helpers/errorHandling'
import {Project, Experiment, Plate, harvestProject} from '../models'
import {QTreeNode} from 'quasar'
import {useI18n} from 'vue-i18n'
import {useRouter} from 'vue-router'
import {useQuasar} from 'quasar'
import {useSettingsStore} from 'stores/settings'
import {storeToRefs} from 'pinia'
import {useProjectStore} from 'stores/project'
import bus from 'src/eventBus'
import {useHarvestStore} from 'stores/harvest-store'

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
  await initialize()
  await harvestStore.initialize()
  bus.on('experiment-updated', () => {
    initialize()
  })
  bus.on('project-updated', () => {
    initialize()
  })
  bus.on('management-command', () => {
    initialize()
  })
})

const harvestStore = useHarvestStore()
const {navigationTree, projectNavigationTree} = storeToRefs(useSettingsStore())
const {harvestProjects} = storeToRefs(harvestStore)
const newProjectDialog = ref<boolean>(false)
const newProjectName = ref<string>('')
const harvestProject = ref<harvestProject | null>(null)
const ownName = ref<boolean>(false)

const options = ref(harvestProjects.value)

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
    // TODO: the handler doesn't work for the 'project' node for some reason
    router.push(`/project/${node.project.id}`)
  } else if ('experiment' in node) {
    router.push(`/project/${node.experiment.project}/experiment/${node.experiment.id}`)
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
  if (project.plates) {
    for (const control_plate of project.plates) {
      addControlPlatesNode(project, control_plate)
    }
  }
  for (const experiment of project.experiments) {
    addExperimentNode(project, experiment)
    for (const plate of experiment.plates) {
      addPlateNode(experiment, plate)
    }
    sortPlateNodes(experiment)
  }
}

const addControlPlatesNode = (project: Project, control_plate: Plate) => {
  const projectNode = projectNodes.value.children?.find(c => c.project.id === project.id)
  if (projectNode) {
    projectNode.children?.push({
      label: `Control plate: ${control_plate.barcode}`,
      icon: 'o_view_module',
      header: 'control_plate',
      handler: nodeHandler,
      plate: control_plate,
    })
  } else {
    handleError(`TSNH: Project ${project.name} not found in tree.`)
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
    // only those children which are experiments and not plates
    const experimentNode = projectNode.children
      ?.filter(c => c.header === 'experiment')
      .find(c => c.experiment.id === experiment.id) // old version:  const experimentNode = projectNode.children?.find(c => c.experiment.id === experiment.id)
    if (experimentNode) {
      experimentNode.children?.push({
        label: `${plate.barcode} (${plate.dimension || t('message.no_dimension')})`,
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
  const projectNode = projectNodes.value.children
    ?.filter(c => c.project !== undefined)
    .find(c => c.project.id === experiment.project)
  if (projectNode) {
    const experimentNode = projectNode.children
      ?.filter(c => c.header === 'experiment')
      .find(c => c.experiment.id === experiment.id)
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
  const payload = {
    name:
      newProjectName.value !== '' && ownName.value
        ? newProjectName.value
        : harvestProject.value
        ? harvestProject.value.name
        : '',
    harvest_id: harvestProject.value ? harvestProject.value.id : null,
    harvest_notes: harvestProject.value ? harvestProject.value.notes : null,
  }
  if (payload.name === '') {
    $q.notify({
      type: 'negative',
      message: t('message.project_name_required'),
    })
    return
  }

  if (projectNodes.value.children) {
    try {
      const project = await projectStore.add(payload)
      addProjectNode(project)
      newProjectName.value = ''
      harvestProject.value = null
    } catch (err) {
      handleError(err)
    }
  }
}

const selectControlLayout = async (project: Project) => {
  await router.push(`/control-plates/${project.id}`)
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

const filterFn = (val: string, update: (arg0: {(): void; (): void}) => void) => {
  if (val === '') {
    update(() => {
      options.value = harvestProjects.value
    })
    return
  }

  update(() => {
    const needle = val.toLowerCase()
    options.value = harvestProjects.value.filter(v => v.name.toLowerCase().indexOf(needle) > -1)
  })
}

const updateHarvestProjects = async () => {
  try {
    await harvestStore.getHarvestProjects()
    await initialize()
  } catch (err) {
    handleError(err)
  }
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
            <q-item-section @click="newProjectDialog = true">{{ t('action.new_project') }}</q-item-section>
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
          <q-separator />
          <q-item clickable v-close-popup>
            <q-item-section @click="selectControlLayout(prop.node.project)">
              {{ t('action.select_control_layout') }}
            </q-item-section>
          </q-item>
        </q-list>
      </q-menu>
      <!--      TODO: node handler doesn't work for 'project' nodes for some reason that is why the @click is here-->
      <span @click="router.push(`/project/${prop.node.project.id}`)">{{ prop.node.label }}</span>
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
  <q-dialog v-model="newProjectDialog" style="width: 700px; max-width: 80vw">
    <q-card>
      <div :style="{display: !ownName && harvestProjects.length > 0 ? 'block' : 'None'}">
        <q-card-section>
          <div class="text-body1">{{ t('message.project_name_harvest') }}:</div>
        </q-card-section>

        <q-card-section class="q-pt-none">
          <q-select
            use-input
            input-debounce="0"
            behavior="dialog"
            v-model="harvestProject"
            :options="options"
            @filter="filterFn"></q-select>

          <q-btn
            class="q-my-lg"
            :label="t('action.update_harvest_projects')"
            icon="update"
            color="secondary"
            @click="updateHarvestProjects" />
        </q-card-section>
      </div>

      <div :style="{display: ownName || harvestProjects.length === 0 ? 'block' : 'None'}">
        <q-card-section>
          <div class="text-body1">{{ t('message.project_name') }}</div>
        </q-card-section>

        <q-card-section class="q-pt-none">
          <q-input
            dense
            v-model="newProjectName"
            autofocus
            :rules="[val => val.length > 0 || 'Please enter something']"></q-input>
        </q-card-section>
      </div>

      <q-toggle
        v-model="ownName"
        :label="t('message.custom_name')"
        right-label
        v-if="harvestProjects.length > 0"></q-toggle>

      <q-card-actions align="right">
        <q-btn flat label="OK" color="primary" v-close-popup @click="newProject"></q-btn>
      </q-card-actions>
    </q-card>
  </q-dialog>
</template>

<style lang="sass">
.q-tree__node-header-content
  cursor: pointer
</style>
