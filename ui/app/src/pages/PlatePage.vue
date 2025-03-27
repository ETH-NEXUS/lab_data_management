<script setup lang="ts">
import {Plate, PlateDimension, Well, WellInfo} from 'src/components/models'
import {ref, onMounted} from 'vue'
import {api} from 'boot/axios'
import {useRoute} from 'vue-router'
import DynamicPlate from 'components/plate/DynamicPlate.vue'
import {handleError} from '../helpers/errorHandling'
import WellDetails from 'components/wells/WellDetails.vue'
import {useI18n} from 'vue-i18n'
import {storeToRefs} from 'pinia'
import {useSettingsStore} from 'stores/settings'
import {PlateLabelValue} from 'components/models'
import {useQuasar} from 'quasar'

const route = useRoute()
const {t} = useI18n()
const $q = useQuasar()

const loading = ref<boolean>(true)
const plate = ref<Plate | null>(null)
const {platePage, projectNavigationTree, libraryNavigationTree, templateNavigationTree} = storeToRefs(
  useSettingsStore()
)
const plateDimensions = ref<Array<PlateDimension>>()

const mapPlateDialog = ref<boolean>(false)
const targetPlateBarcodeOptions = ref<Array<PlateLabelValue>>([])
const applyToAllExperimentPlates = ref<boolean>(false)

const applyTemplateDialog = ref<boolean>(false)
const selectedTemplatePlateId = ref<number>()
const templatePlateBarcodeOptions = ref<Array<PlateLabelValue>>([])
const filteredTemplatePlateOptions = ref<Array<PlateLabelValue>>([])

onMounted(async () => {
  await initialize()
})

const initialize = async () => {
  loading.value = true
  try {
    const resp_plates = await api.get(`/api/plates/?barcode=${route.params.barcode}`)
    if (resp_plates.data.results.length === 1) {
      plate.value = resp_plates.data.results[0]
      platePage.value.selectedWellInfo = undefined
    } else if (resp_plates.data.results.length === 0) {
      handleError(`No plate found with barcode ${route.params.barcode}.`)
    } else {
      handleError(`Multiple plates found with barcode ${route.params.barcode}.`)
    }
    if (!plate.value?.dimension) {
      const resp_dimensions = await api.get('/api/platedimensions/')
      if (resp_dimensions.data.results.length > 0) {
        plateDimensions.value = resp_dimensions.data.results
      }
    }
    const resp_target_plates = await api.get('/api/plates/barcodes/?experiment=true&library=true')
    targetPlateBarcodeOptions.value = resp_target_plates.data.filter(
      (p: PlateLabelValue) => p.label !== plate.value?.barcode
    )
    const resp_template_plates = await api.get(`/api/plates/barcodes/?barcode=${route.params.barcode}`)
    templatePlateBarcodeOptions.value = resp_template_plates.data.filter(
      (p: PlateLabelValue) => p.label !== plate.value?.barcode
    )
  } catch (err) {
    handleError(err)
  } finally {
    loading.value = false
  }
}

const selectedPlateDimension = ref<PlateDimension>()

const setPlateDimension = async () => {
  if (plate.value && selectedPlateDimension.value) {
    try {
      const resp = await api.patch(`/api/plates/${plate.value.id}/`, {
        dimension: selectedPlateDimension.value,
      })
      plate.value = resp.data
      if (plate.value?.experiment) {
        projectNavigationTree.value.needsUpdate = true
      }
      if (plate.value?.library) {
        libraryNavigationTree.value.needsUpdate = true
      }
      if (plate.value?.template) {
        templateNavigationTree.value.needsUpdate = true
      }
    } catch (err) {
      handleError(err)
    }
  } else {
    handleError(t('error.select_plate_dimension'))
  }
}

const wellCreated = (well: Well) => {
  console.log(well)
  console.warn('wellCreated not implemented correctly yet!')
}

const compoundAdded = async (well: Well) => {
  try {
    const resp = await api.get(`/api/wells/${well.id}/`)
    platePage.value.selectedWellInfo = {
      well: resp.data,
      position: resp.data.position,
    }
  } catch (err) {
    handleError(err, false)
  }
}

const measurementAdded = async (well: Well) => {
  try {
    const resp = await api.get(`/api/wells/${well.id}/`)
    platePage.value.selectedWellInfo = {
      well: resp.data,
      position: resp.data.position,
    }
  } catch (err) {
    handleError(err, false)
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
    if (plate.value) {
      $q.loading.show()
      const resp = await api.post(`/api/plates/${plate.value.id}/apply_template/`, {
        template: selectedTemplatePlateId.value,
        apply_to_all_experiment_plates: applyToAllExperimentPlates.value,
      })
      plate.value = resp.data
    }
    $q.loading.hide()
  } catch (err) {
    handleError(err)
  } finally {
    loading.value = false
  }
}

const formatBarcode = (barcode: string | undefined) => {
  if (!barcode) {
    return 'N/A'
  }
  if (barcode.startsWith('__TEMPL__')) {
    return `Templates/${barcode.replace('__TEMPL__', '').replaceAll('_', '/')}`
  }
  return barcode
}

const refresh = async () => {
  await initialize()
}
</script>

<template>
  <q-spinner-grid v-if="loading" color="primary" size="5em" class="absolute-center" />
  <template v-else>
    <q-page v-if="plate?.dimension" class="row items-top q-px-md" :key="`${route.params.barcode}`">
      <q-splitter v-model="platePage.splitter" class="full-width">
        <template v-slot:before>
          <div v-if="plate">
            <h2>{{ formatBarcode(plate.barcode) }}</h2>
            <dynamic-plate
              @refresh="refresh"
              :plate="plate"
              @well-selected="(well_info: WellInfo) => (platePage.selectedWellInfo = well_info)" />
            <div class="q-mt-md q-mb-xl">
              <q-btn
                class="q-mr-xs"
                :label="t('action.map_plate')"
                icon="o_input"
                color="secondary"
                @click="() => (mapPlateDialog = true)" />
              <q-btn
                class="q-ml-md"
                :label="t('action.apply_template')"
                icon="o_input"
                color="secondary"
                @click="() => (applyTemplateDialog = true)" />
            </div>
          </div>
        </template>
        <template v-slot:after>
          <div class="q-px-md" v-if="platePage.selectedWellInfo">
            <well-details
              :key="platePage.selectedWellInfo.position"
              :plate="plate"
              :well-info="platePage.selectedWellInfo"
              @well-created="wellCreated"
              @compound-added="compoundAdded"
              @measurement-added="measurementAdded" />
          </div>
        </template>
      </q-splitter>
    </q-page>
    <div v-else class="fit row wrap justify-left items-start content-start q-ml-md">
      <q-form class="col-6 self-center" @submit="setPlateDimension">
        <h2>{{ formatBarcode(plate?.barcode) }}</h2>
        <q-banner inline-actions rounded class="bg-orange text-white q-mt-lg">
          <span class="text-h6">{{ t('message.plate_has_no_dimension') }}</span>
        </q-banner>
        <q-select
          v-model="selectedPlateDimension"
          :options="plateDimensions"
          option-label="name"
          :label="t('label.plate_dimension')" />
        <div class="fit row justify-end items-start content-start">
          <q-btn
            class="q-mt-xs"
            :disabled="!selectedPlateDimension"
            :label="t('action.submit')"
            type="submit"
            color="primary" />
        </div>
      </q-form>
    </div>
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
          <q-checkbox
            v-model="applyToAllExperimentPlates"
            label="Apply to all experiment plates"></q-checkbox>
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
  </template>
</template>

<style scoped lang="sass">
h2
  font-family: 'Courier New', Courier, monospace
  font-size: 20px

.sticky-header
  /* height or max-height is important */
  height: 310px

  .q-table__top,
  .q-table__bottom,
  thead tr:first-child th
    /* bg color is important for th; just specify one */
    background-color: #c1f4cd

  thead tr th
    position: sticky
    z-index: 1
  thead tr:first-child th
    top: 0

  /* this is when the loading indicator appears */
  &.q-table--loading thead tr:last-child th
    /* height of all previous header rows */
    top: 48px
</style>
