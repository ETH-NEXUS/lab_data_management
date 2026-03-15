<script setup lang="ts">
import { computed, onMounted, ref, watch } from 'vue'
import type { PaginatedResponse } from '~/types/api'
import type { Measurement, MeasurementFeature, Plate, Well, WellInfo } from '~/types/lab'
import { useAPI } from '~/composables/useAPI'
import { hrPositionFromPosition } from '~/utils/plates'
import DynamicImage from '~/components/wells/DynamicImage.vue'
import WellChain from '~/components/wells/WellChain.vue'
import TimeSeriesChart from '~/components/wells/TimeSeriesChartsContainer.vue'
import { usePlateViewStore } from '~/stores/plateView'

const { t } = useI18n()

const addMeasurementDialog = ref<boolean>(false)

const selectedMeasurementFeature = ref<MeasurementFeature>()
const measurementFeatureOptions = ref<MeasurementFeature[]>([])
const filteredMeasurementFeatureOptions = ref<MeasurementFeature[]>([])
const measurementFeatureSearchQuery = ref('')
const enteredMeasurement = ref<number>(0)

const platePage = usePlateViewStore()
const blurCompound = ref<boolean>(false)
const wellInvalid = ref<boolean | undefined>(false)

const props = defineProps<{
  plate: Plate
  wellInfo: WellInfo
}>()

const emit = defineEmits<{
  (e: 'well-created', payload: Well): void
  (e: 'measurement-added', well: WellInfo['well'], measurement: Measurement): void
}>()

const well = ref<Well>()
const loading = ref<boolean>(true)

const isTimeSeries = ref<boolean>(false)

const showTable = computed(() => {
  let result = true
  if (platePage.plotView) {
    result = !isTimeSeries.value
  }
  return result
})

const formatLegacyDateTime = (value: string | Date | undefined): string => {
  if (!value) return 'N/A'
  const date = value instanceof Date ? value : new Date(value)
  if (Number.isNaN(date.getTime())) return String(value)

  const pad = (n: number) => String(n).padStart(2, '0')
  return `${date.getFullYear()}-${pad(date.getMonth() + 1)}-${pad(date.getDate())} ${pad(date.getHours())}:${pad(date.getMinutes())}:${pad(date.getSeconds())}`
}

const initialize = async () => {
  loading.value = true
  try {
    if (props.wellInfo.well) {
      const respWell = await useAPI<Well>(`wells/${props.wellInfo.well.id}/`, { method: 'GET' })
      if (!respWell.error.value && respWell.data.value) {
        well.value = respWell.data.value
      }
    }

    const respMeasurementFeatures = await useAPI<PaginatedResponse<MeasurementFeature>>('measurementfeatures/', {
      method: 'GET',
    })
    if (!respMeasurementFeatures.error.value) {
      measurementFeatureOptions.value = respMeasurementFeatures.data.value?.results ?? []
      filteredMeasurementFeatureOptions.value = measurementFeatureOptions.value
    }
  } catch (err) {
    console.error(err)
  } finally {
    loading.value = false
  }

  wellInvalid.value = well.value?.is_invalid
  if (well.value?.measurements && well.value.measurements.length > 0) {
    const timestamps: Array<string | Date> = []
    well.value.measurements.forEach((measurement) => {
      timestamps.push(measurement.measured_at)
    })
    if ([...new Set(timestamps)].length > 1) {
      isTimeSeries.value = true
    }
  }
}

onMounted(async () => {
  await initialize()
})

watch(
  () => props.wellInfo.well?.id,
  async () => {
    enteredMeasurement.value = 0
    selectedMeasurementFeature.value = undefined
    await initialize()
  },
)

watch(measurementFeatureSearchQuery, (query) => {
  filterMeasurementFeatures(query)
})

const createWell = async () => {
  try {
    const resp = await useAPI<Well>('wells/', {
      method: 'POST',
      body: {
        plate: props.plate.id,
        position: props.wellInfo.position,
      },
    })
    if (!resp.error.value && resp.data.value) {
      emit('well-created', resp.data.value)
    }
  } catch (err) {
    console.error(err)
  }
}

const markWellAsInvalidOrValid = async () => {
  if (!props.wellInfo.well) return
  wellInvalid.value = !wellInvalid.value
  try {
    await useAPI(`wells/${props.wellInfo.well.id}/mark_as_invalid`, { method: 'GET' })
  } catch (err) {
    console.error(err)
  }
}

const addMeasurement = async () => {
  if (selectedMeasurementFeature.value && enteredMeasurement.value !== 0 && props.wellInfo.well) {
    try {
      const resp = await useAPI<Measurement>('measurements/', {
        method: 'POST',
        body: {
          well: props.wellInfo.well.id,
          feature: selectedMeasurementFeature.value.id,
          value: enteredMeasurement.value,
        },
      })
      if (!resp.error.value && resp.data.value) {
        emit('measurement-added', props.wellInfo.well, resp.data.value)
      }
      addMeasurementDialog.value = false
      await initialize()
    } catch (err) {
      console.error(err)
    }
  }
}

const filterMeasurementFeatures = (query: string) => {
  if (query.length > 1) {
    filteredMeasurementFeatureOptions.value = measurementFeatureOptions.value.filter(
      (feature) => feature.name.includes(query) || feature.abbrev.includes(query),
    )
  } else {
    filteredMeasurementFeatureOptions.value = measurementFeatureOptions.value
  }
}

const findAmountFromDonors = () => {
  if (well.value && well.value.donors && well.value.donors.length > 0) {
    let amount = 0
    for (const donor of well.value.donors) {
      amount += donor.amount
    }
    return amount
  }
  return null
}

const onSelectMeasurementFeature = (featureId: string) => {
  selectedMeasurementFeature.value = filteredMeasurementFeatureOptions.value.find(
    (feature) => String(feature.id) === featureId,
  )
}
</script>

<template>
  <template v-if="!loading && well && props.wellInfo.well">
    <div class="container w-full">
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
        <div class="status-box col-4">
          <h4 class="vertical-top">
            {{ t('title.status') }}
          </h4>
        </div>
        <div class="status-box col-8">
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
          <h4 class="vertical-top">
            {{ t('title.amount_dmso') }}
          </h4>

          <div v-if="!wellInvalid" class="text-h6">Valid</div>
          <div v-else class="text-h6">Invalid</div>
          <button type="button" class="action-btn my-sm" @click="markWellAsInvalidOrValid">
            {{ !wellInvalid ? t('action.mark_as_invalid') : t('action.mark_as_valid') }}
          </button>
        </div>
        <div class="col-8">
          <table>
            <tr>
              <th v-if="well.current_info">
                <p>{{ t('well.current_amount') }}: {{ well.current_info.current_amount }} {{ t('unit.mikro') }}</p>
                <p>{{ t('well.current_dmso') }}: {{ well.current_info.current_dmso }}%</p>
              </th>
              <th v-if="findAmountFromDonors()">{{ findAmountFromDonors() }}{{ t('unit.nL') }}</th>
              <th v-if="!well.current_info && !findAmountFromDonors()">
                {{ props.wellInfo.well.amount }}
              </th>
            </tr>
          </table>
        </div>
        <div class="col-12">
          <hr />
        </div>
        <div class="col-12">
          <h4 class="my-none">{{ t('title.compounds') }}</h4>
        </div>
        <div class="mt-md col-12">
          <table>
            <tr>
              <th>{{ t('label.name') }}</th>

              <th>{{ t('label.initial_amount') }}</th>
              <th style="white-space: nowrap">
                {{ t('label.structure') }}
                <input v-model="platePage.showStructure" type="checkbox" @change="blurCompound = !blurCompound" />
              </th>
            </tr>
            <tr v-for="(compound, index) in well.compounds" :key="compound.name + '_' + index">
              <th class="vertical-top">{{ compound.name }}</th>
              <td v-if="compound.amount" class="vertical-top">{{ compound.amount }}{{ t('unit.nL') }}</td>
              <td v-else class="vertical-top">NA</td>
              <td class="vertical-top">
                <DynamicImage
                  v-if="platePage.showStructure && compound.name !== 'unknown'"
                  :key="'' + blurCompound + compound.id"
                  :url="`/api/compounds/${compound.compound}/structure/`"
                  width="250px"
                />
              </td>
            </tr>
          </table>
          <hr />
        </div>
        <div>
          <h4 class="vertical-top">{{ t('title.measurements') }}</h4>
          <div>
            <label v-if="isTimeSeries" class="my-md inline-flex cursor-pointer items-center gap-2 text-sm">
              <input v-model="platePage.plotView" type="checkbox" class="cursor-pointer" />
              <span>{{ t('label.plot_view') }}</span>
            </label>

            <table v-if="showTable" class="mt-lg col-10">
              <thead>
                <tr>
                  <th>{{ t('label.label') }}</th>
                  <th>{{ t('label.timestamp') }}</th>
                  <th>{{ t('label.value') }}</th>
                </tr>
              </thead>
              <tbody>
                <tr v-for="measurement in well.measurements" :key="`${measurement.label}-${measurement.id}`">
                  <td>{{ measurement.label }}</td>
                  <td>{{ measurement.measured_at }}</td>
                  <td>{{ measurement.value }}</td>
                </tr>
              </tbody>
            </table>
          </div>
          <div v-if="well.measurements && well.measurements.length > 0 && platePage.plotView">
            <TimeSeriesChart :measurements="well.measurements" />
          </div>
        </div>
        <div class="col-12">
          <hr />
        </div>
        <div class="col-4">
          <h4 class="vertical-top">
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
              <tr v-for="withdrawal in well.withdrawals" :key="withdrawal.id">
                <td>{{ formatLegacyDateTime(withdrawal.created_at) }}</td>
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
          <h4 class="vertical-top">
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
              <tr v-for="donor in well.donors" :key="donor.id">
                <td>{{ formatLegacyDateTime(donor.created_at) }}</td>
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
          <h4 class="vertical-top">
            {{ t('title.chain') }}
          </h4>
        </div>
        <div class="col-12">
          <WellChain :well="well" />
        </div>
        <div class="col-12">
          <hr />
        </div>
      </div>
    </div>
  </template>
  <template v-if="!props.wellInfo.well">
    <div class="mt-lg bg-orange rounded px-4 py-3 text-white">
      <span class="text-h6">
        {{ t('label.position') }}
        {{ hrPositionFromPosition(props.wellInfo.position, props.plate.dimension) }}:
        {{ t('message.no_well_information') }}
      </span>
      <div class="mt-2">
        <button type="button" class="action-btn action-btn-light" @click="createWell">
          {{ t('action.create_well') }}
        </button>
      </div>
    </div>
  </template>

  <UModal
    :open="addMeasurementDialog"
    :title="t('action.add_measurement')"
    :description="t('hint.measurement_feature')"
    class="w-full sm:max-w-3xl"
    :ui="{ content: 'rounded-2xl bg-white shadow-md' }"
    @update:open="addMeasurementDialog = $event"
  >
    <template #body>
      <div class="space-y-4 p-6">
        <div>
          <label class="mb-1 block text-sm font-medium">{{ t('label.measurement_feature') }}</label>
          <input
            v-model="measurementFeatureSearchQuery"
            type="text"
            class="w-full rounded-full border border-black/15 px-4 py-2 outline-none focus:ring-2 focus:ring-lime-500"
          />
        </div>
        <div>
          <select
            class="w-full rounded-full border border-black/15 px-4 py-2 outline-none focus:ring-2 focus:ring-lime-500"
            @change="onSelectMeasurementFeature(($event.target as HTMLSelectElement).value)"
          >
            <option value="">{{ t('message.no_measurement_feature_found') }}</option>
            <option v-for="feature in filteredMeasurementFeatureOptions" :key="feature.id" :value="feature.id">
              {{ feature.name }} - {{ feature.abbrev }} ({{ feature.unit }})
            </option>
          </select>
        </div>
        <div v-if="selectedMeasurementFeature">
          <label class="mb-1 block text-sm font-medium">{{ t('label.value') }}</label>
          <input
            v-model.number="enteredMeasurement"
            type="number"
            step="0.01"
            class="w-full rounded-full border border-black/15 px-4 py-2 outline-none focus:ring-2 focus:ring-lime-500"
          />
        </div>
      </div>
    </template>
    <template #footer>
      <div class="flex w-full justify-end gap-2 px-6 pb-6">
        <button type="button" class="action-btn action-btn-light" @click="addMeasurementDialog = false">
          {{ t('label.cancel') }}
        </button>
        <button
          type="button"
          class="action-btn"
          :disabled="!selectedMeasurementFeature || enteredMeasurement === 0"
          @click="addMeasurement"
        >
          {{ t('label.add') }}
        </button>
      </div>
    </template>
  </UModal>
</template>

<style scoped>
td,
th {
  text-align: left;
  padding: 5px;
}

h2 {
  font-family: 'Courier New', Courier, monospace;
  font-size: 20px;
}

h4 {
  font-size: 1.5em;
  font-weight: bold;
  margin: 0;
}

.vertical-top {
  vertical-align: top;
}

.text-h6 {
  font-size: 1.25rem;
  font-weight: 500;
}

.bg-orange {
  background: #f97316;
}

.status-box {
  background: #ffedd5;
}

.action-btn {
  cursor: pointer;
  border: 1px solid #134e4a;
  background: transparent;
  color: #134e4a;
  border-radius: 9999px;
  padding: 0.45rem 0.9rem;
  font-size: 0.875rem;
  transition: all 0.2s;
}

.action-btn:hover {
  background: #f0fdfa;
  border-color: #0f766e;
  color: #0f766e;
}

.action-btn:disabled {
  opacity: 0.6;
  cursor: not-allowed;
}

.action-btn-light {
  background: #fff;
  color: #334155;
  border-color: #cbd5e1;
}

.action-btn-light:hover {
  background: #f1f5f9;
  border-color: #94a3b8;
  color: #0f172a;
}

.my-sm {
  margin: 0.35rem 0;
}

.my-none {
  margin: 0;
}

.mt-md {
  margin-top: 1rem;
}

.mt-lg {
  margin-top: 1.25rem;
}

.my-md {
  margin-top: 1rem;
  margin-bottom: 1rem;
}

.col-4 {
  width: 33.333333%;
}

.col-8 {
  width: 66.666667%;
}

.col-10 {
  width: 83.333333%;
}

.col-12 {
  width: 100%;
}

.row {
  display: flex;
  flex-wrap: wrap;
}

.row > .col-4,
.row > .col-8,
.row > .col-10,
.row > .col-12 {
  padding: 0.45rem 0.5rem;
}
</style>
