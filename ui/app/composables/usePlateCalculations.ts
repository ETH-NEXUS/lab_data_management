import { ref } from 'vue'
import {
  PLATE_CALCULATION_DELETE_ERROR_MESSAGE,
  PLATE_CALCULATION_DELETE_PATH,
  PLATE_CALCULATION_ERROR_MESSAGE,
  PLATE_CALCULATION_PATHS,
  PLATE_CALCULATIONS_ENDPOINT,
  PLATE_MEASUREMENT_SOURCES_ERROR_MESSAGE,
  PLATE_MEASUREMENT_SOURCES_PATH,
  type PlateCalculation,
  type PlateCalculationDeleteResponse,
  type PlateCalculationResponse,
  type PlateCalculationSettings,
  type PlateMeasurementSources,
} from '~/types/plateCalculations'
import { requestApiData } from '~/utils/apiRequests'

/**
 * Starts a calculation with the measurements of one plate on the server, or
 * deletes a calculated measurement.
 *
 * Usage examples:
 * - `await calculate(42, 'log10', { label: 'Lum1' })` -> `{ label: 'Lum1_log10', skipped: 0 }`
 * - `await deleteCalculation(42, 'Lum1_log10')` -> `{ label: 'Lum1_log10', deleted: 64 }`
 * - `await fetchSources(42)` -> `{ Lum1: '093026-154654_RKS_300926_1.asc', Lum1_log10: null }`
 */
export const usePlateCalculations = () => {
  const isCalculating = ref(false)

  const calculate = async (
    plateId: number,
    calculation: PlateCalculation,
    settings: PlateCalculationSettings,
  ): Promise<PlateCalculationResponse> => {
    isCalculating.value = true
    try {
      return await requestApiData<PlateCalculationResponse>(
        `${PLATE_CALCULATIONS_ENDPOINT}${plateId}/${PLATE_CALCULATION_PATHS[calculation]}`,
        { method: 'POST', body: settings },
        PLATE_CALCULATION_ERROR_MESSAGE,
      )
    } finally {
      isCalculating.value = false
    }
  }

  const isDeleting = ref(false)

  const deleteCalculation = async (plateId: number, label: string): Promise<PlateCalculationDeleteResponse> => {
    isDeleting.value = true
    try {
      return await requestApiData<PlateCalculationDeleteResponse>(
        `${PLATE_CALCULATIONS_ENDPOINT}${plateId}/${PLATE_CALCULATION_DELETE_PATH}`,
        { method: 'POST', body: { label } },
        PLATE_CALCULATION_DELETE_ERROR_MESSAGE,
      )
    } finally {
      isDeleting.value = false
    }
  }

  const fetchSources = async (plateId: number): Promise<PlateMeasurementSources> => {
    return await requestApiData<PlateMeasurementSources>(
      `${PLATE_CALCULATIONS_ENDPOINT}${plateId}/${PLATE_MEASUREMENT_SOURCES_PATH}`,
      { method: 'GET' },
      PLATE_MEASUREMENT_SOURCES_ERROR_MESSAGE,
    )
  }

  return { isCalculating, calculate, isDeleting, deleteCalculation, fetchSources }
}
