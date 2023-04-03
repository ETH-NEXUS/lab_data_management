<script setup lang="ts">
import {useRoute} from 'vue-router'
import {onMounted, ref} from 'vue'
import {useI18n} from 'vue-i18n'
import {handleError} from 'src/helpers/errorHandling'
import {useProjectStore} from 'stores/project'
import {Project, ProjectPayload} from 'src/components/models'
import {formatDate} from 'src/helpers/dateTime'
import bus from 'src/eventBus'
import {useHarvestStore} from 'stores/harvest-store'

import {useQuasar} from 'quasar'

const route = useRoute()
const projectStore = useProjectStore()
const harvestStore = useHarvestStore()

const $q = useQuasar()
const loading = ref<boolean>(true)
const project = ref<Project | null>(null)

const {t} = useI18n()

const initialize = async () => {
  try {
    await projectStore.initialize()
    const {projects} = projectStore
    project.value = projects?.find((p: Project) => p.id === Number(route.params.project)) ?? null
  } catch (err) {
    handleError(err)
  } finally {
    loading.value = false
  }
}

onMounted(async () => {
  await initialize()
})

const editProject = async (field: string) => {
  if (project.value) {
    $q.dialog({
      title: field === 'name' ? t('project.edit_project_name') : t('project.edit_project_description'),
      message: field === 'name' ? t('project.project_name') : t('project.project_description'),
      prompt: {
        model:
          field === 'name' && project.value.name
            ? project.value.name
            : field === 'description' && project.value.description
            ? project.value.description
            : '',
        type: field === 'name' ? 'text' : 'textarea',
      },
      cancel: true,
      persistent: true,
    }).onOk(async newValue => {
      if (project.value) {
        try {
          const payload = {
            [field]: newValue,
          } as ProjectPayload
          await projectStore.updateProject(project.value.id, payload)
          await initialize()
          bus.emit('project-updated')
        } catch (err) {
          handleError(err, false)
        }
      }
    })
  }
}

const updateHarvestInfo = async () => {
  if (project.value) {
    try {
      const resp = await harvestStore.updateHarvestInfo(project.value.id)
      if (resp.data.success) {
        await initialize()
        $q.notify({
          type: 'positive',
          message: t('message.harvest_info_updated'),
        })
      }
    } catch (err) {
      handleError(err, false)
    }
  }
}
</script>

<template>
  <q-page v-if="project" class="q-px-md">
    <div class="text-h5 q-mt-lg q-mb-md q-pl-xl text-primary">
      {{ t('project.project_name') }}: {{ project.name }}
      <q-btn flat icon="edit" @click="editProject('name')" v-if="!project.harvest_id" />
    </div>

    <div class="q-pa-md row items-start q-gutter-md">
      <q-card class="my-card" flat>
        <q-card-section class="q-pt-xs">
          <div class="text-body1 q-pl-md text-container">
            {{ t('project.project_description') }}:
            <p class="text-body1 text-grey-8">
              {{ project.description || t('project.no_description') }}
              <q-btn flat icon="edit" @click="editProject('description')" />
            </p>
          </div>
          <div class="text-body1 q-pl-md text-container" v-if="project.harvest_notes">
            {{ t('project.harvest_notes') }}:
            <p class="text-body1 text-grey-8">
              {{ project.harvest_notes }}
            </p>
          </div>

          <div class="text-body1 q-pl-md">
            {{ t('project.number_of_experiments') }}:
            <p class="text-body1 text-grey-8">
              {{ project.experiments.length }}
            </p>
          </div>

          <div class="text-body1 q-pl-md">
            {{ t('project.created_at') }}:
            <p class="text-body1 text-grey-8">
              {{ formatDate(project.created_at) }}
            </p>
          </div>

          <q-btn
            v-if="project.harvest_id"
            class="q-ml-md q-mt-md"
            :label="t('message.update_harvest')"
            icon="update"
            color="secondary"
            @click="updateHarvestInfo" />
        </q-card-section>
      </q-card>
    </div>
  </q-page>
</template>

<style scoped lang="sass">

.text-container
  max-width: 600px
  overflow-wrap: anywhere
</style>
