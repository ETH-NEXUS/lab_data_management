import { ref } from 'vue'
import {
  PLATE_CALCULATION_ERROR_MESSAGE,
  PLATE_CALCULATION_PATHS,
  PLATE_CALCULATIONS_ENDPOINT,
  type PlateCalculation,
  type PlateCalculationResponse,
  type PlateCalculationSettings,
} from '~/types/plateCalculations'
import { requestApiData } from '~/utils/apiRequests'

/**
 * Starts a calculation with the measurements of one plate on the server.
 *
 * Usage example:
 * - `await calculate(42, 'log10', { label: 'Lum1' })` -> `{ label: 'Lum1_log10', skipped: 0 }`
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

  return { isCalculating, calculate }
}
