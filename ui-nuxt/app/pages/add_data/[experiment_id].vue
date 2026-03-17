<script setup lang="ts">
import { computed, watch } from 'vue'
import { useExperimentStore } from '~/stores/experiments'
import { getErrorMessage } from '~/utils/errors'

const route = useRoute()
const toast = useToast()
const { t } = useI18n()
const experimentStore = useExperimentStore()

/**
 * Current experiment id parsed from route param.
 *
 * Route param examples:
 * - `/add_data/21` -> `21`
 * - `/add_data/not-a-number` -> `null`
 */
const experimentId = computed<number | null>(() => {
  const rawParam = route.params.experiment_id
  const value = Array.isArray(rawParam) ? rawParam[0] : rawParam
  const parsedId = Number(value)

  if (!Number.isInteger(parsedId) || parsedId <= 0) {
    return null
  }

  return parsedId
})

const hasValidExperimentId = computed<boolean>(() => experimentId.value !== null)

/**
 * Loads prefilled plate info rows for the current experiment.
 *
 * Loaded data example:
 * - `[{ measurement_label: 'OD600', measurement_timestamp: ['2026-01-01T00:00:00Z'], plate_barcode: 'A001', lib_plate_barcode: 'LIB-1', replicate: '1', cell_type: 'HeLa', condition: 'Control' }]`
 */
const loadPrefilledPlateInfo = async () => {
  if (!experimentId.value) return
  await experimentStore.fetchPrefilledPlateInfo(experimentId.value)
}

watch(
  experimentId,
  () => {
    void loadPrefilledPlateInfo()
  },
  { immediate: true },
)

/**
 * Downloads current table rows as CSV file.
 *
 * CSV header example:
 * - `measurement_label,measurement_timestamp,plate,lib_plate_barcode,replicate,cell_type,condition`
 */
const downloadCsv = () => {
  const columns = [
    'measurement_label',
    'measurement_timestamp',
    'plate_barcode',
    'lib_plate_barcode',
    'replicate',
    'cell_type',
    'condition',
  ] as const

  let headerLine = columns.join(',') + '\n'
  headerLine = headerLine.replace('plate_barcode', 'plate')

  const rowLines: string[] = []

  for (const row of experimentStore.prefilledPlateInfo) {
    const cells: string[] = []

    for (const column of columns) {
      const value = row[column]

      if (Array.isArray(value)) {
        cells.push(`"${value.join(',')}"`)
      } else {
        cells.push(`"${value}"`)
      }
    }

    rowLines.push(cells.join(','))
  }

  const csv = headerLine + rowLines.join('\n')
  const blob = new Blob([csv], { type: 'text/csv' })
  const url = window.URL.createObjectURL(blob)
  const link = document.createElement('a')

  link.href = url
  link.download = 'experiment_data.csv'
  link.click()

  window.URL.revokeObjectURL(url)
}

const savePlateInfo = async () => {
  if (!experimentId.value) return

  try {
    await experimentStore.savePlateInfo({
      experiment_id: experimentId.value,
      plate_info: experimentStore.prefilledPlateInfo,
    })

    toast.add({
      title: t('experiments.add_data.save_success'),
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: t('experiments.add_data.save_failed'),
      description: getErrorMessage(err),
      color: 'error',
      duration: 3500,
    })
  }
}
</script>

<template>
  <section class="space-y-5 p-6">
    <p v-if="!hasValidExperimentId" class="text-sm text-slate-600">
      {{ t('experiments.add_data.invalid_experiment_id') }}
    </p>

    <template v-else>
      <h1 class="text-primary text-center text-2xl font-semibold">
        {{ t('experiments.add_data.title') }}
      </h1>

      <UCard
        class="my-card mx-auto"
        :ui="{
          root: 'border border-white/40 bg-white/30 backdrop-blur-md shadow-sm divide-y divide-white/20',
        }"
      >
        <p v-if="experimentStore.isLoadingPrefilledPlateInfo" class="text-sm text-slate-600">
          {{ t('experiments.add_data.loading') }}
        </p>

        <p v-else-if="experimentStore.prefilledPlateInfo.length === 0" class="text-sm text-slate-600">
          {{ t('experiments.add_data.empty') }}
        </p>

        <ClientOnly v-else>
          <div class="text-body1">
            <vue-excel-editor
              v-model="experimentStore.prefilledPlateInfo"
              :page="20"
              :no-paging="false"
              :filter-row="true"
              :autocomplete="true"
              :autocomplete-count="100"
              :readonly="false"
              :width="'100%'"
              :spellcheck="true"
              :new-if-bottom="true"
              :remember="true"
              :enter-to-south="true"
              :disable-panel-setting="false"
              :disable-panel-filter="false"
              :no-mouse-scroll="false"
            >
              <vue-excel-column
                :readonly="true"
                field="measurement_label"
                label="Measurement Label"
                type="string"
                width="150px"
              />
              <vue-excel-column
                :readonly="true"
                field="measurement_timestamp"
                label="Measurement Timestamps"
                type="string"
                width="200px"
              />
              <vue-excel-column
                :readonly="true"
                field="plate_barcode"
                label="Plate Barcode"
                type="string"
                width="120px"
              />
              <vue-excel-column
                :readonly="false"
                field="lib_plate_barcode"
                label="Library Plate Barcode"
                type="string"
                width="150px"
              />
              <vue-excel-column field="replicate" label="Replicate" type="string" width="100px" />
              <vue-excel-column field="cell_type" label="Cell Type" type="string" width="100px" />
              <vue-excel-column field="condition" label="Condition" type="string" width="100px" />
            </vue-excel-editor>
          </div>
        </ClientOnly>

        <template #footer>
          <div class="flex flex-wrap gap-2">
            <UButton
              color="primary"
              :label="t('experiments.add_data.save_button')"
              :loading="experimentStore.isSavingPlateInfo"
              @click="savePlateInfo"
            />
            <UButton
              color="neutral"
              variant="outline"
              :label="t('experiments.add_data.download_csv')"
              @click="downloadCsv"
            />
          </div>
        </template>
      </UCard>
    </template>
  </section>
</template>

<style scoped>
.my-card {
  max-width: 1100px;
}
</style>
