<script setup lang="ts">
import {useRoute} from 'vue-router'
import {onMounted, ref, computed} from 'vue'
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
  PlateLabelValue,
} from 'src/components/models'
import {formatDate} from 'src/helpers/dateTime'
import GenerateBarcodeForm from '../components/GenerateBarcodeForm.vue'
import {downloadCSVData, generateBarcodes} from 'components/helpers'
import {csvColumnsNames} from 'components/data'
import {useQuasar} from 'quasar'
import {api} from 'boot/axios'
import bus from 'src/eventBus'
import MeasurementCalculator from 'components/MeasurementCalculator.vue'
import ExperimentHeatmap from 'components/ExperimentHeatmap.vue'
import {storeToRefs} from 'pinia'
import {useSettingsStore} from 'stores/settings'

const route = useRoute()
const projectStore = useProjectStore()

const $q = useQuasar()
const options = ref<DimensionsOption[]>([])
const dimension = ref<number | null>(null)
const loading = ref<boolean>(true)
const project = ref<Project | null>(null)
const experiment = ref<Experiment | null>(null)
const generateBarcodeDialogToggle = ref<boolean>(false)
const addNewMeasurementDialog = ref<boolean>(false)

const applyTemplateDialog = ref<boolean>(false)
const selectedTemplatePlateId = ref<number>()
const templatePlateBarcodeOptions = ref<Array<PlateLabelValue>>([])
const filteredTemplatePlateOptions = ref<Array<PlateLabelValue>>([])

const {t} = useI18n()
const {showExperimentResults} = storeToRefs(useSettingsStore())

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
  const resp_template_plates = await api.get('/api/plates/barcodes/?template=true')
  templatePlateBarcodeOptions.value = resp_template_plates.data
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
        message: t('message.select_dimension'),
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

const filterTemplatePlates = (query: string, update: (f: () => void) => void) => {
  update(() => {
    if (query.length > 1) {
      filteredTemplatePlateOptions.value = templatePlateBarcodeOptions.value.filter(m =>
        m.label.includes(query)
      )
    }
    filteredTemplatePlateOptions.value = templatePlateBarcodeOptions.value.sort((a, b) =>
      a.label.localeCompare(b.label)
    )
  })
}

const applyTemplate = async () => {
  try {
    if (experiment.value) {
      $q.loading.show({
        message: t('info.applying_in_progress'),
      })
      await api.post(`/api/experiments/bulk_apply_template/`, {
        template: selectedTemplatePlateId.value,
        experiment_id: experiment.value.id,
      })
      $q.loading.hide()
    }
  } catch (err) {
    handleError(err)
  } finally {
    loading.value = false
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
</script>

<template>
  <q-page v-if="project" class="q-px-md">
    <div v-if="experiment" class="text-h5 q-mt-lg q-mb-md q-pl-xl text-primary">
      {{ project.name }}: {{ getExperiment(Number(route.params.experiment)).name }}
      <q-btn flat icon="edit" @click="editExperiment('name')" />
    </div>
    <div class="q-pa-md row items-start q-gutter-md">
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

          <div class="text-body1 text-grey-8" v-if="experiment.barcode_specifications">
            <q-list>
              <div v-for="(s, i) in experiment.barcode_specifications" :key="`${s.prefix}-${s.id}`">
                <q-item>
                  <q-item-section>
                    <!--                    <q-item-label class="text-body1 q-mt-lg">Barcode set #{{ i + 1 }}</q-item-label>-->
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
                          :label="t('action.add_plates')"
                          type="submit"
                          color="secondary"
                          class="q-mt-md"
                          @click="addPlatesToExperiment(experiment.id, s.id)"></q-btn>
                      </div>
                      <div v-else class="text-subtitle2 text-red-3">* {{ t('experiment.plates_added') }}</div>
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
                  :label="t('action.delete_specification')"
                  color="red"
                  class="q-mt-md"
                  @click="deleteBarcode(s.id)" />

                <q-btn
                  flat
                  :label="t('action.edit_specification')"
                  color="warning"
                  class="q-mt-md"
                  @click="openEditField(i)" />

                <q-btn
                  flat
                  v-if="experiment.barcode_specifications.length > 0"
                  :label="t('action.download_csv')"
                  type="submit"
                  color="secondary"
                  class="q-mt-md"
                  @click="download" />

                <div :id="`edit-${i}`" :class="`hidden q-mt-lg q-ml-md`">
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

        <q-card-actions class="q-ml-md">
          <q-btn
            class="q-ml-xs"
            :label="t('action.add_barcode_specification')"
            icon="qr_code_2"
            color="secondary"
            @click="generateBarcodeDialogToggle = true" />
          <q-btn
            class="q-ml-xs"
            :label="t('action.apply_template')"
            icon="o_layers"
            color="secondary"
            @click="applyTemplateDialog = true" />
          <q-btn
            v-if="experiment.available_measurement_labels.length > 0"
            class="q-ml-xs"
            :label="t('action.add_measurement')"
            icon="o_layers"
            color="secondary"
            @click="addNewMeasurement" />
        </q-card-actions>

        <q-card-section class="q-mt-lg">
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
              :available-measurement-labels="experiment.available_measurement_labels"
              :overall-stats="experiment.details.overall_stats"
              :experiment-id="experiment.id" />
          </q-expansion-item>
        </q-card-section>
      </q-card>
    </div>
    <q-dialog v-model="generateBarcodeDialogToggle">
      <GenerateBarcodeForm :experiment-id="experiment.id" @update="update" />
    </q-dialog>
    <q-dialog v-model="applyTemplateDialog" persistent>
      <q-card style="width: 700px; max-width: 80vw" class="q-px-sm">
        <q-card-body class="q-gutter-y-sm">
          <q-select
            filled
            v-model="selectedTemplatePlateId"
            emit-value
            map-options
            use-input
            input-debounce="0"
            :label="t('label.template_plate')"
            :options="filteredTemplatePlateOptions"
            @filter="filterTemplatePlates"
            behavior="menu"
            :hint="t('hint.template_plate')">
            <template v-slot:no-option>
              <q-item>
                <q-item-section class="text-grey">
                  {{ t('message.no_plates_found') }}
                </q-item-section>
              </q-item>
            </template>
          </q-select>
        </q-card-body>
        <q-card-actions align="right" class="bg-white text-teal">
          <q-btn flat :label="t('label.cancel')" v-close-popup />
          <q-btn
            flat
            :label="t('action.apply')"
            :disabled="!selectedTemplatePlateId"
            v-close-popup
            @click="applyTemplate" />
        </q-card-actions>
      </q-card>
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
  </q-page>
</template>

<style scoped lang="sass">

.hidden
  visibility: hidden

.calculator_dialog
  min-width: 800px
</style>
