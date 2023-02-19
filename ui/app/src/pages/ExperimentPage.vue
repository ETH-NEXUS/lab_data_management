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
const editToggle = ref<boolean>(false)

const {t} = useI18n()

const initialize = async () => {
  try {
    await projectStore.initialize()
    if (projectStore.projects) {
      project.value =
        projectStore.projects.find((p: Project) => p.id === Number(route.params.project)) || null
      experiment.value = await getExperiment(Number(route.params.experiment))
    }
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
const deleteBarcode = async (id: number) => {
  try {
    await projectStore.deleteBarcode(id)
    await initialize()
  } catch (err) {
    handleError(err)
  }
}

const openEditField = (index: number) => {
  const id = `edit-${index}`
  const editForm = document.getElementById(id) as HTMLInputElement
  if (editForm) {
    if (editForm.classList.contains('hidden')) {
      editForm.classList.remove('hidden')
    } else {
      editForm.classList.add('hidden')
    }
  }
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
          <div class="text-overline">{{t('description')}}:</div>
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

          <div class="text-overline">Barcode sets:</div>

          <div
            class="text-body1 text-grey-8"
            v-if="getExperiment(Number(route.params.experiment)).barcode_specifications">
            <q-list>
              <div v-for="(s, i) in experiment.barcode_specifications" :key="`${s.prefix}-${s.id}`">
                <q-item>
                  <q-item-section>
                    <q-item-label class="text-subtitle1 q-mt-lg">Barcode set #{{ i + 1 }}</q-item-label>

                    <q-item-label caption lines="2">
                      <q-table
                        :rows="generateBarcodes(s.prefix, s.number_of_plates, s.sides)"
                        row-key="name"></q-table>
                    </q-item-label>
                  </q-item-section>
                </q-item>
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

                <q-btn
                  flat
                  label="Delete specification"
                  color="red"
                  class="q-mt-md"
                  @click="deleteBarcode(s.id)"></q-btn>

                <q-btn
                  flat
                  label="Edit specification"
                  color="warning"
                  class="q-mt-md"
                  @click="openEditField(i)"></q-btn>

                <q-separator class="q-my-md"></q-separator>
                <div>>> Add paltes to the experiment using these specifications/div>

                <div :id="`edit-${i}`" :class="`hidden q-mt-lg`">
                  <GenerateBarcodeForm
                    :edit="true"
                    @update="update"
                    :experiment-id="experiment.id"
                    :prefilledData="{
                      index: i,
                      prefix: s.prefix,
                      number_of_plates: s.number_of_plates,
                      sides: s.sides,
                      id: s.id,
                    }" />
                </div>
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

<style scoped lang="sass">

.hidden
  visibility: hidden
</style>
