import { defineStore } from 'pinia'
import { ref } from 'vue'
import type { PlateWellContentMode, PlatePaletteOption } from '~/types/plates'
import type { WellInfo } from '~/types/lab'
import { getDefaultPlatePalette } from '~/utils/plateHeatmap'

export const usePlateViewStore = defineStore('plateViewStore', () => {
  const splitter = ref(58)
  const selectedWellInfo = ref<WellInfo | undefined>(undefined)
  const measurementOptions = ref<string[] | null>(null)
  const selectedMeasurement = ref<string | null>(null)
  const selectedTimestampIdx = ref(0)
  const selectedPosControl = ref<string | null>(null)
  const selectedNegControl = ref<string | null>(null)
  const wellContent = ref<PlateWellContentMode>('hr_position')
  const showHeatmap = ref(false)
  const smallerMapView = ref(false)
  const heatmapPalette = ref<PlatePaletteOption>(getDefaultPlatePalette())
  const plotView = ref(false)
  const showStructure = ref(true)
  const squareCompoundType = ref(false)

  /**
   * Resets only view-local selection state for a newly opened plate.
   *
   * Reset examples:
   * - `selectedWellInfo = undefined`
   * - `showHeatmap = false`
   * - `wellContent = 'hr_position'`
   */
  const resetForNewPlate = (): void => {
    selectedWellInfo.value = undefined
    measurementOptions.value = null
    selectedMeasurement.value = null
    selectedTimestampIdx.value = 0
    selectedPosControl.value = null
    selectedNegControl.value = null
    wellContent.value = 'hr_position'
    showHeatmap.value = false
    smallerMapView.value = false
    heatmapPalette.value = getDefaultPlatePalette()
    plotView.value = false
    showStructure.value = true
    squareCompoundType.value = false
  }

  return {
    splitter,
    selectedWellInfo,
    measurementOptions,
    selectedMeasurement,
    selectedTimestampIdx,
    selectedPosControl,
    selectedNegControl,
    wellContent,
    showHeatmap,
    smallerMapView,
    heatmapPalette,
    plotView,
    showStructure,
    squareCompoundType,
    resetForNewPlate,
  }
})
