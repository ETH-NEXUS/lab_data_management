<script setup lang="ts">
import {Plate, PlateDimension, Well, WellInfo, Compound} from 'src/components/models'
import {ref, onMounted} from 'vue'
import {api} from '../boot/axios'
import {useRoute} from 'vue-router'
import DynamicPlate from '../components/DynamicPlate.vue'
import {handleError} from '../helpers/errorHandling'
import WellDetails from '../components/WellDetails.vue'
import {useI18n} from 'vue-i18n'
import {storeToRefs} from 'pinia'
import {useSettingsStore} from '../stores/settings'

const route = useRoute()
const {t} = useI18n()

const loading = ref<boolean>(true)
const plate = ref<Plate | null>(null)
const {platePage} = storeToRefs(useSettingsStore())
const plateDimensions = ref<Array<PlateDimension>>()

onMounted(async () => {
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
  } catch (err) {
    handleError(err)
  } finally {
    loading.value = false
  }
})

const selectedPlateDimension = ref<PlateDimension>()

const setPlateDimension = async () => {
  if (plate.value && selectedPlateDimension.value) {
    try {
      const resp = await api.patch(`/api/plates/${plate.value.id}/`, {
        dimension: selectedPlateDimension.value,
      })
      plate.value = resp.data
    } catch (err) {
      handleError(err)
    }
  } else {
    handleError(t('error.select_plate_dimension'))
  }
}

const wellCreated = (well: Well) => {
  plate.value?.wells?.push(well)
  platePage.value.selectedWellInfo = {
    well: well,
    position: well.position,
  }
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
</script>

<template>
  <q-spinner-grid v-if="loading" color="primary" size="5em" class="absolute-center" />
  <template v-else>
    <q-page v-if="plate?.dimension" class="row items-top q-px-md" :key="`${route.params.barcode}`">
      <q-splitter v-model="platePage.splitter" class="full-width">
        <template v-slot:before>
          <div v-if="plate">
            <h2>{{ plate.barcode }}</h2>
            <dynamic-plate
              :plate="plate"
              @well-selected="(well_info: WellInfo) => (platePage.selectedWellInfo = well_info)" />
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
        <h2>{{ plate?.barcode }}</h2>
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
  </template>
</template>

<style scoped lang="sass">
h2
  font-family: 'Courier New', Courier, monospace
</style>
