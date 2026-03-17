<script setup lang="ts">
import { computed, ref } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import ExperimentBarcodeSpecificationTable from '~/components/experiments/ExperimentBarcodeSpecificationTable.vue'
import { useExperimentStore } from '~/stores/experiments'
import { BARCODE_CSV_COLUMNS, downloadCsvData, generateBarcodes } from '~/utils/barcode'
import type { Barcode, BarcodeSpecification } from '~/types/lab'
import type { Experiment } from '~/types/experiments'
import { getErrorMessage } from '~/utils/errors'

const props = defineProps<{
  experiment: Experiment
}>()

const emit = defineEmits<{
  (e: 'add-specification'): void
  (e: 'deleted-specification', specificationId: number): void
}>()

const { t } = useI18n()
const toast = useToast()
const experimentStore = useExperimentStore()
const deletingSpecificationId = ref<number | null>(null)

/**
 * Normalized barcode specification list for rendering.
 *
 * Returned data example:
 * - `[{ id: 1, prefix: 'A001', number_of_plates: 12, sides: ['North', 'South'], experiment: 21, created_at: '...', modified_at: '...' }]`
 */
const barcodeSpecifications = computed<BarcodeSpecification[]>(() => props.experiment.barcode_specifications ?? [])

/**
 * Flattens generated barcodes from all specifications into one CSV export list.
 *
 * Returned data example:
 * - `[{ NorthBarcode: 'A001_1$A001_1', SouthBarcode: '', EastBarcode: '', WestBarcode: '' }]`
 */
const allGeneratedBarcodes = computed<Barcode[]>(() => {
  return barcodeSpecifications.value.flatMap((specification) =>
    generateBarcodes(specification.prefix, specification.number_of_plates, specification.sides),
  )
})

const canDownloadBarcodesCsv = computed(() => allGeneratedBarcodes.value.length > 0)

const onAddSpecification = () => emit('add-specification')

const formatSidesLabel = (sides: string[]): string => {
  if (sides.length === 0) return '-'
  return sides.join(', ')
}

/**
 * Downloads one combined CSV for all generated barcode rows.
 *
 * Output file example:
 * - `barcodes.csv`
 */
const onDownloadBarcodesCsv = () => {
  if (!canDownloadBarcodesCsv.value) return

  downloadCsvData(BARCODE_CSV_COLUMNS, allGeneratedBarcodes.value, 'barcodes.csv')
}

/**
 * Deletes one barcode specification and notifies the parent page to refresh queries.
 *
 * Accepted id example:
 * - `14`
 */
const onDeleteSpecification = async (specificationId: number) => {
  deletingSpecificationId.value = specificationId

  try {
    await experimentStore.deleteBarcodeSpecification(specificationId)
    emit('deleted-specification', specificationId)

    toast.add({
      title: t('experiments.barcodes.delete_success'),
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: t('experiments.barcodes.delete_failed'),
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  } finally {
    deletingSpecificationId.value = null
  }
}
</script>

<template>
  <UCard
    class="mx-auto w-[80%]"
    :ui="{
      root: 'core-card divide-y divide-slate-200/70',
    }"
  >
    <template #header>
      <div class="flex flex-wrap items-center justify-between gap-2">
        <p class="font-semibold">{{ t('experiments.barcodes.title') }}</p>

        <div class="flex flex-wrap items-center gap-2">
          <BaseButton
            v-if="canDownloadBarcodesCsv"
            :label="t('experiments.barcodes.download_button')"
            :on-click="onDownloadBarcodesCsv"
            variant="outline"
            size="sm"
            width="auto"
            class-name="text-blue-700 hover:text-blue-800 hover:bg-blue-50 hover:border-blue-200"
          />

          <BaseButton
            :label="t('experiments.barcodes.add_button')"
            :on-click="onAddSpecification"
            variant="outline"
            size="sm"
            width="auto"
            class-name="text-blue-700 hover:text-blue-800 hover:bg-blue-50 hover:border-blue-200"
          />
        </div>
      </div>
    </template>

    <div v-if="barcodeSpecifications.length === 0" class="text-sm text-slate-600">
      {{ t('experiments.barcodes.empty') }}
    </div>

    <div v-else class="space-y-3">
      <div v-for="specification in barcodeSpecifications" :key="specification.id" class="space-y-2">
        <details class="group rounded-xl border border-slate-200 bg-white/90 p-3">
          <summary class="flex cursor-pointer list-none items-center justify-between gap-3">
            <div class="min-w-0 text-sm text-slate-700">
              <p class="font-semibold text-slate-800">{{ specification.prefix }}</p>
              <p class="truncate text-xs text-slate-600">
                {{ t('experiments.barcodes.fields.number_of_plates') }}: {{ specification.number_of_plates }} |
                {{ t('experiments.barcodes.fields.sides') }}: {{ formatSidesLabel(specification.sides) }}
              </p>
            </div>

            <UIcon
              name="i-heroicons-chevron-right"
              class="size-5 shrink-0 text-slate-500 transition-transform duration-200 group-open:rotate-90"
            />
          </summary>

          <div class="mt-3 space-y-2">
            <ExperimentBarcodeSpecificationTable :specification="specification" />

            <BaseButton
              :label="t('experiments.barcodes.delete_button')"
              :on-click="() => onDeleteSpecification(specification.id)"
              variant="outline"
              size="sm"
              width="auto"
              class-name="text-red-700 hover:text-red-800 hover:bg-red-50 hover:border-red-200"
              :loading="deletingSpecificationId === specification.id"
              :disabled="deletingSpecificationId !== null && deletingSpecificationId !== specification.id"
            />
          </div>
        </details>
      </div>
    </div>
  </UCard>
</template>
