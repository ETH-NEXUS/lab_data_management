<script setup lang="ts">
import { computed, watch } from 'vue'
import DynamicPlate from '~/components/plates/DynamicPlate.vue'
import PlatePaneHeader from '~/components/plates/PlatePaneHeader.vue'
import ResizableSplitPane from '~/components/common/ResizableSplitPane.vue'
import WellDetails from '~/components/wells/WellDetails.vue'
import {
  PLATE_PAGE_DEFAULT_LEFT_PERCENT,
  PLATE_PAGE_MAX_LEFT_PERCENT,
  PLATE_PAGE_MIN_LEFT_PERCENT,
} from '~/types/plates'
import { usePlateStore } from '~/stores/plates'
import { usePlateViewStore } from '~/stores/plateView'
import { formatPlateBarcodeLabel, getPlateRouteIdLabel } from '~/utils/plates'
import { useAPI } from '~/composables/useAPI'
import type { Well, WellDetails as WellDetailsType, WellInfo } from '~/types/lab'

const route = useRoute()
const { t } = useI18n()
const plateStore = usePlateStore()
const plateViewStore = usePlateViewStore()

/**
 * Current barcode route param used for API plate lookup.
 *
 * Returned data examples:
 * - `'demo_1'`
 * - `'250401control'`
 * - `'-'` (missing route fallback)
 */
const plateRouteBarcode = computed(() => getPlateRouteIdLabel(route.params.id))

const splitPercent = computed({
  get: () => plateViewStore.splitter || PLATE_PAGE_DEFAULT_LEFT_PERCENT,
  set: (value: number) => {
    plateViewStore.splitter = value
  },
})

const pageTitle = computed(() => {
  const barcode = plateStore.currentPlate?.barcode ?? plateRouteBarcode.value
  return t('plates.page.title', { id: formatPlateBarcodeLabel(barcode) })
})

const loadPlatePage = async (): Promise<void> => {
  const barcode = plateRouteBarcode.value
  if (!barcode || barcode === '-') return

  plateViewStore.resetForNewPlate()
  await plateStore.initializePlatePage(barcode)
}

watch(
  () => route.fullPath,
  () => {
    void loadPlatePage()
  },
  { immediate: true },
)

const onWellSelected = (wellInfo: WellInfo): void => {
  plateViewStore.selectedWellInfo = wellInfo
}

const wellCreated = async (_well: Well) => {
  await loadPlatePage()
}

const refreshSelectedWell = async (wellId: number | undefined) => {
  if (!wellId) return
  try {
    const resp = await useAPI<Well | WellDetailsType>(`wells/${wellId}/`, { method: 'GET' })
    if (!resp.error.value && resp.data.value) {
      const refreshedWell = resp.data.value
      plateViewStore.selectedWellInfo = {
        well: refreshedWell as unknown as WellDetailsType,
        position: refreshedWell.position,
      }
    }
  } catch (err) {
    console.error(err)
  }
}

const measurementAdded = async (well: WellInfo['well']) => {
  await refreshSelectedWell(well?.id)
}
</script>

<template>
  <section class="mb-50 h-dvh overflow-y-auto pb-100">
    <ResizableSplitPane
      v-model="splitPercent"
      :min-left-percent="PLATE_PAGE_MIN_LEFT_PERCENT"
      :max-left-percent="PLATE_PAGE_MAX_LEFT_PERCENT"
      :divider-width="3"
    >
      <template #left>
        <div class="h-full p-4">
          <div class="h-full">
            <PlatePaneHeader
              :title="pageTitle"
              :description="t('plates.page.left_description')"
              title-class="mb-3 font-mono text-3xl font-medium text-slate-800"
            />

            <div v-if="plateStore.isInitializingPlatePage" class="flex min-h-[320px] items-center justify-center">
              <div
                class="flex flex-col items-center gap-3 rounded-2xl border border-teal-100 bg-white/70 px-8 py-6 shadow-sm backdrop-blur-sm"
              >
                <span
                  class="inline-block h-12 w-12 animate-spin rounded-full border-[3px] border-teal-900/25 border-t-teal-900"
                />
                <p class="text-sm font-medium text-teal-900">
                  {{ t('plates.page.loading') }}
                </p>
              </div>
            </div>
            <p v-else-if="plateStore.error" class="rounded-lg border border-red-200 bg-red-50 p-3 text-sm text-red-700">
              {{ plateStore.error }}
            </p>
            <p v-else-if="!plateStore.currentPlate" class="text-sm text-slate-600">
              {{ t('plates.page.not_found') }}
            </p>
            <DynamicPlate
              v-else
              :plate="plateStore.currentPlate"
              @well-selected="onWellSelected"
              @refresh="loadPlatePage"
            />
          </div>
        </div>
      </template>

      <template #right>
        <div class="h-full p-4">
          <div class="h-full">
            <PlatePaneHeader :title="t('plates.page.right_title')" :description="t('plates.page.right_description')" />

            <div v-if="plateStore.currentPlate && plateViewStore.selectedWellInfo" class="mt-4">
              <WellDetails
                :key="plateViewStore.selectedWellInfo.position"
                :plate="plateStore.currentPlate"
                :well-info="plateViewStore.selectedWellInfo"
                @well-created="wellCreated"
                @measurement-added="measurementAdded"
              />
            </div>
            <p v-else class="mt-4 text-sm text-slate-600">
              {{ t('plates.page.right_placeholder') }}
            </p>
          </div>
        </div>
      </template>
    </ResizableSplitPane>
  </section>
</template>
