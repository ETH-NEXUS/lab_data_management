<script setup lang="ts">
import {useRoute} from 'vue-router'
import {onMounted, ref} from 'vue'
import {useI18n} from 'vue-i18n'
import {handleError} from 'src/helpers/errorHandling'
import {useProjectStore} from 'stores/project'
import {Project, ProjectPayload} from 'src/components/models'
import {formatDate} from 'src/helpers/dateTime'

import {useQuasar} from 'quasar'

const route = useRoute()
const projectStore = useProjectStore()

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
  $q.dialog({
    title: field === 'name' ? t('project.edit_project_name') : t('project.edit_project_description'),
    message: field === 'name' ? t('project.project_name') : t('project.project_description'),
    prompt: {
      model: '',
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
      } catch (err) {
        handleError(err, false)
      }
    }
  })
}
</script>

<template>
  <q-page v-if="project" class="q-px-md">
    <div class="text-h5 q-mt-lg q-mb-md q-pl-xl text-primary">
      {{ t('project.project_name') }}: {{ project.name }}
      <q-btn flat icon="edit" @click="editProject((field = 'name'))" />
    </div>
    <div class="q-pa-md row items-start q-gutter-md">
      <q-card class="my-card" flat>
        <q-card-section class="q-pt-xs">
          <div class="text-body1 q-pl-md">
            {{ t('project.project_description') }}:
            <p class="text-body1 text-grey-8">
              {{ project.description || t('project.no_description') }}
              <q-btn flat icon="edit" @click="editProject((field = 'description'))" />
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
        </q-card-section>
      </q-card>
    </div>
  </q-page>
</template>
