<script setup lang="ts">
import {useRoute, useRouter} from 'vue-router'
import {onMounted, ref} from 'vue'
import {useI18n} from 'vue-i18n'
import {handleError} from 'src/helpers/errorHandling'
import {useProjectStore} from 'stores/project'
import {
  Barcode,
  DimensionsOption,
  Experiment,
  PlateDimension,
  Project,
  ExperimentPayload,
} from 'src/components/models'
import {formatDate} from 'src/helpers/dateTime'
import GenerateBarcodeForm from 'components/experiment/GenerateBarcodeForm.vue'
import {downloadCSVData, generateBarcodes} from 'components/helpers'
import {csvColumnsNames} from 'components/data'
import {useQuasar} from 'quasar'
import bus from 'src/eventBus'
import MeasurementCalculator from 'components/plate/MeasurementCalculator.vue'
import ExperimentHeatmap from 'components/experiment/ExperimentHeatmap.vue'
import {storeToRefs} from 'pinia'
import {useSettingsStore} from 'stores/settings'
import GenerateReportDialog from 'components/experiment/GenerateReportDialog.vue'
import DownloadCSVDialog from 'components/experiment/DownloadCSVDialog.vue'

const route = useRoute()
const router = useRouter()
const projectStore = useProjectStore()

const $q = useQuasar()
const options = ref<DimensionsOption[]>([])
const loading = ref<boolean>(true)
const project = ref<Project | null>(null)
const experiment = ref<Experiment | null>(null)
const generateBarcodeDialogToggle = ref<boolean>(false)
const addNewMeasurementDialog = ref<boolean>(false)
const expanded = ref<boolean>(false)
const generateReportDialog = ref<boolean>(false)
const downloadCSVDialog = ref<boolean>(false)

const {t} = useI18n()
const {showExperimentResults} = storeToRefs(useSettingsStore())

const initialize = async () => {
  try {
    await projectStore.initialize()
    const {projects, plateDimensions, experiments} = projectStore
    project.value = projects?.find((p: Project) => p.id === Number(route.params.project)) ?? null
    options.value = plateDimensions?.map((d: PlateDimension) => ({label: d.name, value: d.id})) ?? []
    experiment.value = experiments?.find((e: Experiment) => e.id === Number(route.params.experiment)) ?? null
    showExperimentResults.value = false
    if (experiment.value) {
      await projectStore.getNotebookOutputFiles(experiment.value.name)
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
  generateBarcodeDialogToggle.value = false
}
const deleteBarcode = async (id: number) => {
  try {
    await projectStore.deleteBarcode(id)
    await initialize()
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

const editExperiment = async (field: string) => {
  if (experiment.value) {
    $q.dialog({
      title: field === 'name' ? t('title.edit_experiment_name') : t('title.edit_experiment_description'),
      message: field === 'name' ? t('message.experiment_name') : t('message.experiment_description'),
      prompt: {
        model:
          field === 'name' && experiment.value.name
            ? experiment.value.name
            : field === 'description' && experiment.value.description
            ? experiment.value.description
            : '',
        type: field === 'name' ? 'text' : 'textarea',
      },
      cancel: true,
      persistent: true,
    }).onOk(async newValue => {
      if (experiment.value) {
        try {
          const payload = {
            [field]: newValue,
          } as ExperimentPayload
          await projectStore.updateExperiment(experiment.value.id, payload)
          await initialize()
          bus.emit('experiment-updated')
        } catch (err) {
          handleError(err, false)
        }
      }
    })
  }
}

const calculateNewMeasurement = async (expression: string, newLabel: string, usedLabels: string[]) => {
  addNewMeasurementDialog.value = false
  $q.loading.show({
    message: t('info.calculation_in_progress'),
  })
  await projectStore.addNewMeasurement(null, expression, newLabel, usedLabels, experiment.value?.id)
  $q.loading.hide()
}

const addNewMeasurement = () => {
  addNewMeasurementDialog.value = true
}

const generateReport = async () => {
  generateReportDialog.value = true
}

const downloadReport = async (path: string) => {
  await projectStore.downloadPDFReport(path)
}

const downloadCsvData = async () => {
  downloadCSVDialog.value = true
}
</script>

<template>
  <q-page v-if="project" class="q-px-md">
    <div v-if="experiment" class="text-h5 q-mt-lg q-mb-md q-pl-xl text-primary">
      {{ project.name }}: {{ getExperiment(Number(route.params.experiment)).name }}
      <q-btn flat icon="edit" @click="editExperiment('name')" />
    </div>

    <div class="q-pa-md row items-start q-gutter-md" v-if="experiment">
      <q-card class="my-card" flat>
        <q-card-section class="q-pt-xs">
          <div class="text-body1 q-pl-md">
            {{ t('experiment.description') }}:
            <p class="text-body1 text-grey-8">
              {{ experiment.description || t('info.no_description') }}
              <q-btn flat icon="edit" @click="editExperiment('description')" />
            </p>
          </div>

          <div class="text-body1 q-pl-md">
            {{ t('experiment.number_plates') }}:
            <p class="text-body1 text-grey-8">
              {{ experiment.plates.length }}
            </p>
          </div>

          <div class="text-body1 q-pl-md">
            {{ t('experiment.created_at') }}:
            <p class="text-body1 text-grey-8">
              {{ formatDate(experiment.created_at) }}
            </p>
          </div>

          <div class="text-body1 q-pl-md">{{ t('experiment.barcode_sets') }}:</div>

          <div class="text-body1 text-grey-8">
            <q-list>
              <div v-for="(s, i) in experiment.barcode_specifications" :key="`${s.prefix}-${s.id}-${i}`">
                <q-item>
                  <q-item-section>
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
                  :label="t('action.delete_specification')"
                  color="red"
                  class="q-mt-md"
                  @click="deleteBarcode(s.id)" />

                <q-btn
                  flat
                  v-if="experiment.barcode_specifications.length > 0"
                  :label="t('action.download_csv')"
                  type="submit"
                  color="secondary"
                  class="q-mt-md"
                  @click="download" />
              </div>
            </q-list>
          </div>
          <q-btn flat color="primary" @click="router.push(`/add_data/${experiment.id}`)">
            Add Experiment Data
          </q-btn>
        </q-card-section>

        <q-card-actions class="q-ml-md">
          <q-btn
            class="q-ml-xs"
            :label="t('action.add_barcode_specification')"
            icon="qr_code_2"
            color="secondary"
            @click="generateBarcodeDialogToggle = true" />

          <q-btn
            v-if="experiment.available_measurement_labels.length > 0"
            class="q-ml-xs"
            :label="t('action.add_measurement')"
            icon="o_layers"
            color="secondary"
            @click="addNewMeasurement" />
          <q-btn
            v-if="experiment.available_measurement_labels.length > 0"
            class="q-ml-xs"
            :label="t('action.generate_report')"
            icon="o_layers"
            color="secondary"
            @click="generateReport" />
          <q-btn
            v-if="experiment.available_measurement_labels.length > 0"
            class="q-ml-xs"
            :label="t('action.download_csv_data')"
            icon="o_layers"
            color="secondary"
            @click="downloadCsvData" />
        </q-card-actions>
        <!--  @click="generateReport" />-->
        <q-card-section v-if="projectStore.outputNotebooks.length > 0" class="q-ml-md">
          <p class="text-caption subtitle">Available reports:</p>
          <p
            @click="downloadReport(notebook)"
            class="text-blue-8 cursor-pointer"
            :key="notebook"
            v-for="notebook in projectStore.outputNotebooks">
            {{ notebook.split('/').slice(-1)[0] }}
          </p>
        </q-card-section>

        <q-card-section
          class="q-mt-md"
          v-if="experiment.details.measurement_labels && experiment.details.measurement_labels.length > 0">
          <q-expansion-item
            v-model="expanded"
            class="shadow-1 overflow-hidden"
            style="border-radius: 30px"
            icon="explore"
            :label="t('action.show_results')"
            @show="showExperimentResults = true"
            @hide="showExperimentResults = false"
            header-class="bg-secondary text-white"
            expand-icon-class="text-white">
            <ExperimentHeatmap
              v-if="showExperimentResults"
              :timestamps="experiment.details.measurement_timestamps"
              :available-measurement-labels="experiment.details.measurement_labels"
              :overall-stats="experiment.details.overall_stats"
              :experiment-id="experiment.id" />
          </q-expansion-item>
        </q-card-section>
      </q-card>
    </div>
    <div v-if="experiment">
      <q-dialog v-model="generateBarcodeDialogToggle">
        <GenerateBarcodeForm :experiment-id="experiment.id" @update="update" />
      </q-dialog>
      <q-dialog v-model="addNewMeasurementDialog">
        <q-card class="calculator_dialog">
          <q-card-section>
            <MeasurementCalculator
              :labels="experiment.available_measurement_labels"
              @calculate="calculateNewMeasurement" />
          </q-card-section>
        </q-card>
      </q-dialog>
      <q-dialog v-model="generateReportDialog">
        <GenerateReportDialog
          :labels="experiment.available_measurement_labels"
          :label="experiment.available_measurement_labels[0]"
          :experiment-name="experiment.name" />
      </q-dialog>
      <q-dialog v-model="downloadCSVDialog">
        <DownloadCSVDialog
          :labels="experiment.available_measurement_labels"
          :experiment-name="experiment.name"
          :label="experiment.available_measurement_labels[0]" />
      </q-dialog>
    </div>
  </q-page>
</template>

<style scoped lang="sass">

.hidden
  visibility: hidden

.calculator_dialog
  min-width: 800px

.subtitle
  font-size: 1.1rem
</style>
