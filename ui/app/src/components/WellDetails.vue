<script setup lang="ts">
import {defineProps, defineEmits, PropType, onMounted} from 'vue'
import {Plate, WellInfo, Well, MeasurementFeature, Compound} from './models'
import DynamicImage from './DynamicImage.vue'
import {useI18n} from 'vue-i18n'
import {api} from '../boot/axios'
import {handleError} from '../helpers/errorHandling'
import {hrPositionFromPosition} from '../helpers/plate'
import {ref} from 'vue'
import {date} from 'quasar'
import {storeToRefs} from 'pinia'
import {useSettingsStore} from '../stores/settings'
import WellChain from './WellChain.vue'

const {t} = useI18n()

const addCompoundDialog = ref<boolean>(false)
const addMeasurementDialog = ref<boolean>(false)

const selectedCompound = ref<Compound>()
const compoundOptions = ref<Array<Compound>>([])

const selectedWell = ref<Well>()
const wellOptions = ref<Array<Well>>([])
const filteredWellOptions = ref<Array<Well>>([])
const selectedAmount = ref<number>(0)

const selectedMeasurementFeature = ref<MeasurementFeature>()
const measurementFeatureOptions = ref<Array<MeasurementFeature>>([])
const filteredMeasurementFeatureOptions = ref<Array<MeasurementFeature>>([])
const enteredMeasurement = ref<number>(0)

const {wellDetails} = storeToRefs(useSettingsStore())
const blurCompound = ref<boolean>(false)

const props = defineProps({
  plate: {
    type: Object as PropType<Plate>,
    required: true,
  },
  wellInfo: {
    type: Object as PropType<WellInfo>,
    required: true,
  },
})

const emit = defineEmits(['well-created', 'compound-added', 'measurement-added'])

onMounted(async () => {
  try {
    const resp = await api.get('/api/measurementfeatures/')
    measurementFeatureOptions.value = resp.data.results
  } catch (err) {
    handleError(err, false)
  }
})

const createWell = async () => {
  try {
    const resp = await api.post('/api/wells/', {
      plate: props.plate.id,
      position: props.wellInfo.position,
    })
    emit('well-created', resp.data as Well)
  } catch (err) {
    handleError(err, false)
  }
}

const searchCompounds = (query: string, update: (f: () => void) => void) => {
  if (query.length > 2) {
    update(async () => {
      try {
        const resp = await api.get(`/api/compounds/?search=${query}`)
        compoundOptions.value = resp.data.results
      } catch (err) {
        handleError(err, false)
      }
    })
  }
}

const addCompound = async () => {
  if (selectedCompound.value && selectedWell.value && selectedAmount.value > 0) {
    try {
      await api.post('/api/wellcompounds/', {
        well: props.wellInfo.well.id,
        compound: selectedCompound.value.id,
        amount: selectedAmount.value,
      })
      await api.post('/api/wellwithdrawals/', {
        well: selectedWell.value.id,
        target_well: props.wellInfo.well.id,
        amount: selectedAmount.value,
      })
      emit('compound-added', props.wellInfo.well, selectedCompound.value)
    } catch (err) {
      handleError(err, false)
    }
  }
}

const addMeasurement = async () => {
  if (selectedMeasurementFeature.value && enteredMeasurement.value != 0) {
    try {
      const resp = await api.post('/api/measurements/', {
        well: props.wellInfo.well.id,
        feature: selectedMeasurementFeature.value.id,
        value: enteredMeasurement.value,
      })
      const measurement = resp.data
      emit('measurement-added', props.wellInfo.well, measurement)
    } catch (err) {
      handleError(err, false)
    }
  }
}

const updateWellOptions = () => {
  if (selectedCompound.value) {
    wellOptions.value = []
    for (const well of selectedCompound.value.wells) {
      wellOptions.value.push(well)
    }
  }
}

const filterWells = (query: string, update: (f: () => void) => void) => {
  update(() => {
    if (query.length > 1) {
      filteredWellOptions.value = wellOptions.value.filter(
        w => w.plate.library !== null && !w.mixture && w.plate.barcode.includes(query)
      )
    } else {
      filteredWellOptions.value = wellOptions.value.filter(w => w.plate.library !== null && !w.mixture)
    }
  })
}

const filterMeasurementFeatures = (query: string, update: (f: () => void) => void) => {
  update(() => {
    if (query.length > 1) {
      filteredMeasurementFeatureOptions.value = measurementFeatureOptions.value.filter(
        m => m.name.includes(query) || m.abbrev.includes(query)
      )
    } else {
      filteredMeasurementFeatureOptions.value = measurementFeatureOptions.value
    }
  })
}
</script>

<template>
  <template v-if="props.wellInfo?.well">
    <div class="container full-width">
      <div class="row">
        <div class="col-12">
          <h2>
            {{ props.wellInfo.well.hr_position }} ({{ props.wellInfo.well.type }})
            <small style="font-size: 16px">({{ props.wellInfo.well.position }})</small>
            <small style="font-size: 10px">[{{ props.wellInfo.well.id }}]</small>
          </h2>
        </div>
      </div>
      <div v-if="props.wellInfo.well.status" class="row">
        <div class="col-4 bg-orange-2">
          <h4 class="q-ma-none vertical-top">
            <q-icon name="o_warning" color="red" />
            {{ t('title.status') }}
          </h4>
        </div>
        <div class="col-8 bg-orange-2">
          <table>
            <tr>
              <th>{{ props.wellInfo.well.status }}</th>
            </tr>
          </table>
        </div>
        <div class="col-12">
          <hr />
        </div>
      </div>
      <div class="row">
        <div class="col-4">
          <h4 class="q-ma-none vertical-top">
            <q-icon name="o_water_drop" />
            {{ t('title.amount') }}
          </h4>
        </div>
        <div class="col-8">
          <table>
            <tr>
              <th>{{ props.wellInfo.well.amount }}{{ t('unit.amount') }}</th>
            </tr>
          </table>
        </div>
        <div class="col-12">
          <hr />
        </div>
        <div class="col-12">
          <h4 class="q-my-none">{{ t('title.compounds') }}</h4>
          <q-btn
            :label="t('action.add_compound')"
            icon="o_colorize"
            color="secondary"
            @click="addCompoundDialog = true" />
        </div>
        <div class="col-12 q-mt-md">
          <table>
            <tr>
              <th>{{ t('label.name') }}</th>
              <th>{{ t('label.identifier') }}</th>
              <th>{{ t('label.initial_amount') }}</th>
              <th style="white-space: nowrap">
                {{ t('label.structure') }}
                <q-toggle
                  v-model="wellDetails.showStructure"
                  @update:model-value="blurCompound = !blurCompound" />
              </th>
            </tr>
            <tr v-for="compound in props.wellInfo.well.compounds" :key="compound.identifier">
              <th class="vertical-top">{{ compound.name }}</th>
              <td class="vertical-top">{{ compound.identifier }}</td>
              <td class="vertical-top">{{ compound.amount }}{{ t('unit.amount') }}</td>
              <td class="vertical-top">
                <dynamic-image
                  v-if="wellDetails.showStructure"
                  :url="`/api/compounds/${compound.compound}/structure/`"
                  width="250px"
                  :key="'' + blurCompound + compound.id" />
              </td>
            </tr>
          </table>
          <hr />
        </div>
        <div class="col-4">
          <h4 class="q-ma-none vertical-top">{{ t('title.measurements') }}</h4>
          <q-btn
            :label="t('action.add_measurement')"
            icon="o_straighten"
            color="secondary"
            @click="addMeasurementDialog = true" />
        </div>
        <div class="col-8">
          <table>
            <thead>
              <tr>
                <th>{{ t('label.name') }}</th>
                <th>{{ t('label.abbrev') }}</th>
                <th>{{ t('label.value') }}</th>
                <th>{{ t('label.unit') }}</th>
              </tr>
            </thead>
            <tbody>
              <tr v-for="measurement in props.wellInfo.well.measurements" :key="measurement.feature.name">
                <td>{{ measurement.feature.name }}</td>
                <td>{{ measurement.feature.abbrev }}</td>
                <td>{{ measurement.value }}</td>
                <td>{{ measurement.feature.unit }}</td>
              </tr>
            </tbody>
          </table>
        </div>
        <div class="col-12">
          <hr />
        </div>
        <div class="col-4">
          <h4 class="q-ma-none vertical-top">
            <q-icon name="remove" />
            {{ t('title.withdrawals') }}
          </h4>
        </div>
        <div class="col-8">
          <table>
            <thead>
              <tr>
                <th>{{ t('label.timestamp') }}</th>
                <th>{{ t('label.amount') }}</th>
                <th>{{ t('label.target_well') }}</th>
              </tr>
            </thead>
            <tbody>
              <tr v-for="withdrawal in props.wellInfo.well.withdrawals" :key="withdrawal.id">
                <td>{{ date.formatDate(withdrawal.created_at, 'YYYY-MM-DD hh:mm:ss') }}</td>
                <td>{{ withdrawal.amount }}{{ t('unit.amount') }}</td>
                <td v-if="withdrawal.target_well">
                  {{ withdrawal.target_well.plate.barcode }}: {{ withdrawal.target_well.hr_position }}
                </td>
                <td v-else>N/A</td>
              </tr>
            </tbody>
          </table>
        </div>
        <div class="col-12">
          <hr />
        </div>
        <div class="col-4">
          <h4 class="q-ma-none vertical-top">
            <q-icon name="add" />
            {{ t('title.donors') }}
          </h4>
        </div>
        <div class="col-8">
          <table>
            <thead>
              <tr>
                <th>{{ t('label.timestamp') }}</th>
                <th>{{ t('label.amount') }}</th>
                <th>{{ t('label.source_well') }}</th>
              </tr>
            </thead>
            <tbody>
              <tr v-for="donor in props.wellInfo.well.donors" :key="donor.id">
                <td>{{ date.formatDate(donor.created_at, 'YYYY-MM-DD hh:mm:ss') }}</td>
                <td>{{ donor.amount }}{{ t('unit.amount') }}</td>
                <td v-if="donor.well">{{ donor.well.plate.barcode }}: {{ donor.well.hr_position }}</td>
                <td v-else>N/A</td>
              </tr>
            </tbody>
          </table>
        </div>
        <div class="col-12">
          <hr />
        </div>
        <div class="col-12">
          <h4 class="q-ma-none vertical-top">
            <q-icon name="o_route" />
            {{ t('title.chain') }}
          </h4>
        </div>
        <div class="col-12">
          <well-chain :well="props.wellInfo.well" />
        </div>
        <div class="col-12">
          <hr />
        </div>
      </div>
    </div>
  </template>
  <template v-else>
    <q-banner inline-actions rounded class="bg-orange text-white q-mt-lg">
      <span class="text-h6">
        {{ t('label.position') }}
        {{ hrPositionFromPosition(props.wellInfo.position, props.plate.dimension) }}:
        {{ t('message.no_well_information') }}
      </span>

      <template v-slot:action>
        <q-btn flat :label="t('action.create_well')" @click="createWell" />
      </template>
    </q-banner>
  </template>
  <q-dialog v-model="addCompoundDialog" persistent>
    <q-card style="width: 700px; max-width: 80vw" class="q-px-sm">
      <q-card-body class="q-gutter-y-sm">
        <q-select
          filled
          v-model="selectedCompound"
          use-input
          input-debounce="0"
          :label="t('label.compound')"
          :options="compoundOptions"
          option-label="name"
          @filter="searchCompounds"
          @update:model-value="updateWellOptions"
          behavior="menu"
          :hint="t('hint.compound_to_add')">
          <template v-slot:no-option>
            <q-item>
              <q-item-section class="text-grey">
                {{ t('message.no_compounds_found') }}
              </q-item-section>
            </q-item>
          </template>
        </q-select>
        <q-select
          v-if="selectedCompound"
          filled
          v-model="selectedWell"
          use-input
          input-debounce="0"
          :label="t('label.well')"
          :options="filteredWellOptions"
          :option-label="
            w => `${w.plate.library} - ${w.plate.barcode}: ${w.hr_position} (${w.amount}${t('unit.amount')})`
          "
          @filter="filterWells"
          behavior="menu"
          :hint="t('hint.well_to_transfer_from')">
          <template v-slot:no-option>
            <q-item>
              <q-item-section class="text-grey">
                {{ t('message.no_wells_found') }}
              </q-item-section>
            </q-item>
          </template>
        </q-select>
        <q-input
          v-if="selectedWell"
          filled
          v-model="selectedAmount"
          :label="t('label.amount')"
          mask="#.##"
          fill-mask="0"
          reverse-fill-mask
          :hint="t('hint.amount_to_transfer')"
          input-class="text-right" />
      </q-card-body>
      <q-card-actions align="right" class="bg-white text-teal">
        <q-btn flat :label="t('label.cancel')" v-close-popup />
        <q-btn
          flat
          :label="t('label.add')"
          :disabled="!selectedCompound || !selectedWell || selectedAmount == 0"
          v-close-popup
          @click="addCompound" />
      </q-card-actions>
    </q-card>
  </q-dialog>
  <q-dialog v-model="addMeasurementDialog" persistent>
    <q-card style="width: 700px; max-width: 80vw" class="q-px-sm">
      <q-card-body class="q-gutter-y-sm">
        <q-select
          filled
          v-model="selectedMeasurementFeature"
          use-input
          input-debounce="0"
          :label="t('label.measurement_feature')"
          :options="filteredMeasurementFeatureOptions"
          :option-label="m => `${m.name} - ${m.abbrev} (${m.unit})`"
          @filter="filterMeasurementFeatures"
          behavior="menu"
          :hint="t('hint.measurement_feature')">
          <template v-slot:no-option>
            <q-item>
              <q-item-section class="text-grey">
                {{ t('message.no_measurement_feature_found') }}
              </q-item-section>
            </q-item>
          </template>
        </q-select>
        <q-input
          v-if="selectedMeasurementFeature"
          filled
          v-model="enteredMeasurement"
          :label="t('label.value')"
          mask="#.##"
          fill-mask="0"
          reverse-fill-mask
          :hint="t('hint.measurement_value')"
          input-class="text-right" />
      </q-card-body>
      <q-card-actions align="right" class="bg-white text-teal">
        <q-btn flat :label="t('label.cancel')" v-close-popup />
        <q-btn
          flat
          :label="t('label.add')"
          :disabled="!selectedMeasurementFeature || enteredMeasurement == 0"
          v-close-popup
          @click="addMeasurement" />
      </q-card-actions>
    </q-card>
  </q-dialog>
</template>

<style scoped lang="sass">
td, th
  text-align: left
h2
  font-family: 'Courier New', Courier, monospace
  font-size: 20px
h4
  font-size: 1.5em
  font-weight: bold
</style>
