<script setup lang="ts">
import {useRoute} from 'vue-router'
import {onMounted, ref} from 'vue'
import {useI18n} from 'vue-i18n'
import {handleError} from 'src/helpers/errorHandling'
import {useProjectStore} from 'stores/project'
import {Project, Experiment} from 'src/components/models'
import {formatDate} from 'src/helpers/dateTime'
import GenerateBarcodeForm from '../components/GenerateBarcodeForm.vue'
import {downloadCSVData, generateBarcodes} from 'components/helpers'
import {csvColumnsNames} from 'components/data'

const route = useRoute()
const projectStore = useProjectStore()

const loading = ref<boolean>(true)
const project = ref<Project | null>(null)
const experiment = ref<Experiment | null>(null)
const generateBarcodeDialogToggle = ref<boolean>(false)

const {t} = useI18n()

const initialize = async () => {
  try {
    await projectStore.initialize()
    if (projectStore.projects) {
      project.value =
        projectStore.projects.find((p: Project) => p.id === Number(route.params.project)) || null
      experiment.value = await getExperiment(Number(route.params.experiment))
    }
    console.log('projectStore', projectStore.projects)
  } catch (err) {
    handleError(err)
  } finally {
    loading.value = false
  }
}

onMounted(async () => {
  await initialize()
})

const getExperiment = (id: number) => {
  return project.value?.experiments.find((e: Experiment) => e.id === id) || null
}

const update = async () => {
  await initialize()
}
</script>

<template>
  <q-page v-if="project" class="q-px-md">
    <div v-if="experiment" class="text-h5 q-mt-lg q-mb-md text-primary">
      {{ project.name }}: {{ getExperiment(Number(route.params.experiment)).name }}
    </div>
    <div class="q-pa-md row items-start q-gutter-md">
      <q-card class="my-card" flat>
        <q-card-section class="q-pt-xs">
          <div class="text-overline">Description:</div>
          <div class="text-body1 text-grey-8">
            {{ experiment.description || 'No description provided' }}
          </div>
          <div class="text-overline">Number of plates:</div>
          <div class="text-body1 text-grey-8">
            {{ experiment.plates.length }}
          </div>
          <div class="text-overline">Created at:</div>
          <div class="text-body1 text-grey-8">
            {{ formatDate(experiment.created_at) }}
          </div>

          <div
            class="text-body1 text-grey-8 q-mt-lg"
            v-if="getExperiment(Number(route.params.experiment)).barcode_specifications">
            <q-list>
              <div v-for="(s, i) in experiment.barcode_specifications" :key="`${s.prefix}-${s.id}`">
                <q-item>
                  <q-item-section>
                    <q-item-label class="text-secondary">Barcode set #{{ i + 1 }}</q-item-label>
                    <q-btn
                      flat
                      label="Download csv"
                      type="submit"
                      color="secondary"
                      class="q-mt-md"
                      @click="
                        downloadCSVData(
                          csvColumnsNames,
                          generateBarcodes(s.prefix, s.number_of_plates, s.sides),
                          'barcodes.csv'
                        )
                      "></q-btn>
                    <q-item-label caption lines="2">
                      <q-table
                        :rows="generateBarcodes(s.prefix, s.number_of_plates, s.sides)"
                        row-key="name"></q-table>
                    </q-item-label>
                  </q-item-section>
                </q-item>
              </div>
            </q-list>
          </div>
        </q-card-section>

        <q-card-actions>
          <q-btn color="secondary" class="q-mt-lg" @click="generateBarcodeDialogToggle = true">
            Add barcodes
          </q-btn>
        </q-card-actions>
      </q-card>
    </div>
    <q-dialog v-model="generateBarcodeDialogToggle">
      <GenerateBarcodeForm :experiment-id="experiment.id" @update="update" />
    </q-dialog>
  </q-page>
</template>

<style scoped lang="sass"></style>
