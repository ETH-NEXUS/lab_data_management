<script setup lang="ts">
import {useRoute} from 'vue-router'
import {onMounted, ref} from 'vue'
import {useI18n} from 'vue-i18n'
import {handleError} from 'src/helpers/errorHandling'
import {useProjectStore} from 'stores/project'
import {Project, Experiment, DimensionsOption, Barcode} from 'src/components/models'
import {formatDate} from 'src/helpers/dateTime'
import GenerateBarcodeForm from '../components/GenerateBarcodeForm.vue'
import {downloadCSVData, generateBarcodes} from 'components/helpers'
import {csvColumnsNames} from 'components/data'
import {PlateDimension} from 'src/components/models'
import {useQuasar} from 'quasar'

const route = useRoute()
const projectStore = useProjectStore()

const $q = useQuasar()
const options = ref<DimensionsOption[]>([])
const dimension = ref<number | null>(null)
const loading = ref<boolean>(true)
const project = ref<Project | null>(null)
const experiment = ref<Experiment | null>(null)
const generateBarcodeDialogToggle = ref<boolean>(false)

const {t} = useI18n()

const initialize = async () => {
  try {
    await projectStore.initialize()
    const {projects, plateDimensions, experiments} = projectStore
    project.value = projects?.find((p: Project) => p.id === Number(route.params.project)) ?? null
    options.value = plateDimensions?.map((d: PlateDimension) => ({label: d.name, value: d.id})) ?? []
    experiment.value = experiments?.find((e: Experiment) => e.id === Number(route.params.experiment)) ?? null
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
    if (editForm) {
      editForm.classList.toggle('hidden')
    }
  }
}

const openDimensionsOptions = (index: number) => {
  const id = `dimensions-${index}`
  const dimensionsForm = document.getElementById(id) as HTMLInputElement
  if (dimensionsForm) {
    dimensionsForm.classList.toggle('hidden')
  }
}

const addPlatesToExperiment = async (experimentId: number, barcodeSpecificationsId: number) => {
  try {
    if (dimension.value) {
      await projectStore.addPlatesToExperiment(experimentId, barcodeSpecificationsId, dimension.value)
      await initialize()
    } else {
      $q.notify({
        type: 'negative',
        message: 'Please select a dimension',
      })
    }
  } catch (err) {
    handleError(err)
  }
}

const download = (): void => {
  let barcodes: Barcode[] = []
  const barcodeSpecifications = experiment.value?.barcode_specifications || []

  for (const spec of barcodeSpecifications) {
    const generatedBarcodes = generateBarcodes(spec.prefix, spec.number_of_plates, spec.sides)
    barcodes = [...barcodes, ...generatedBarcodes]
  }

  downloadCSVData(csvColumnsNames, barcodes, 'barcodes.csv')
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
          <div class="text-overline">{{ t('experiment.description') }}:</div>
          <div class="text-body1 text-grey-8">
            {{ experiment.description || 'No description provided' }}
          </div>
          <div class="text-overline">{{ t('experiment.number_plates') }}:</div>
          <div class="text-body1 text-grey-8">
            {{ experiment.plates.length }}
          </div>
          <div class="text-overline">{{ t('experiment.created_at') }}:</div>
          <div class="text-body1 text-grey-8">
            {{ formatDate(experiment.created_at) }}
          </div>

          <div class="text-overline">{{ t('experiment.barcode_sets') }}:</div>

          <q-btn
            label="Download csv"
            type="submit"
            color="secondary"
            class="q-mt-md"
            @click="download"></q-btn>

          <div
            class="text-body1 text-grey-8"
            v-if="getExperiment(Number(route.params.experiment)).barcode_specifications">
            <q-list>
              <div v-for="(s, i) in experiment.barcode_specifications" :key="`${s.prefix}-${s.id}`">
                <q-item>
                  <q-item-section>
                    <q-item-label class="text-body1 q-mt-lg">Barcode set #{{ i + 1 }}</q-item-label>
                    <div
                      class="text-primary cursor-pointer q-mt-lg q-mb-sm"
                      @click="openDimensionsOptions(i)">
                      >> {{ t('experiment.add_plates') }}
                    </div>
                    {{}}
                    <div :id="`dimensions-${i}`" :class="`hidden q-my-lg`">
                      <div
                        v-if="
                          !(
                            experiment.plates.length > 0 &&
                            experiment.plates[0].barcode &&
                            experiment.plates[0].barcode.includes(s.prefix)
                          )
                        ">
                        <p class="text-subtitle2">{{ t('experiment.choose_dimensions') }}:</p>
                        <q-option-group :options="options" type="radio" v-model="dimension"></q-option-group>

                        <q-btn
                          label="Add plates"
                          type="submit"
                          color="secondary"
                          class="q-mt-md"
                          @click="addPlatesToExperiment(experiment.id, s.id)"></q-btn>
                      </div>
                      <div v-else class="text-subtitle2 text-red-3">
                        *You have already added the plates with these barcode specifications to this
                        experiment
                      </div>
                    </div>

                    <q-item-label caption lines="2">
                      <q-table
                        bordered
                        :rows="generateBarcodes(s.prefix, s.number_of_plates, s.sides)"
                        row-key="name"></q-table>
                    </q-item-label>
                  </q-item-section>
                </q-item>

                <q-btn
                  flat
                  label="Delete specifications"
                  color="red"
                  class="q-mt-md"
                  @click="deleteBarcode(s.id)"></q-btn>

                <q-btn
                  flat
                  label="Edit specifications"
                  color="warning"
                  class="q-mt-md"
                  @click="openEditField(i)"></q-btn>

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
            Add barcode specifications
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
