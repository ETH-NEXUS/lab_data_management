import { ref } from 'vue'
import {
  BACKGROUND_CORRECTION_ENDPOINT,
  LOG10_ENDPOINT,
  LOG10_ERROR_MESSAGE,
  type Log10Response,
} from '~/types/backgroundCorrection'
import { requestApiData } from '~/utils/apiRequests'

/**
 * Starts the log10 of a measurement of one plate on the server.
 *
 * Usage example:
 * - `await log10Measurement(42, 'Lum1')` -> `{ label: 'Lum1_log10', skipped: 2 }`
 */
export const usePlateLog10 = () => {
  const isCalculatingLog10 = ref(false)

  const log10Measurement = async (plateId: number, label: string): Promise<Log10Response> => {
    isCalculatingLog10.value = true
    try {
      return await requestApiData<Log10Response>(
        `${BACKGROUND_CORRECTION_ENDPOINT}${plateId}/${LOG10_ENDPOINT}`,
        { method: 'POST', body: { label } },
        LOG10_ERROR_MESSAGE,
      )
    } finally {
      isCalculatingLog10.value = false
    }
  }

  return { isCalculatingLog10, log10Measurement }
}
