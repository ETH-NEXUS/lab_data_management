import { ref } from 'vue'
import {
  ACTIVITY_ENDPOINT,
  ACTIVITY_ERROR_MESSAGE,
  BACKGROUND_CORRECTION_ENDPOINT,
  type ActivitySettings,
  type BackgroundCorrectionResponse,
} from '~/types/backgroundCorrection'
import { requestApiData } from '~/utils/apiRequests'

/**
 * Starts the %Activity of a measurement of one plate on the server.
 *
 * Usage example:
 * - `await calculateActivity(42, { label: 'Lum1', negative_type: 'N', positive_type: 'P' })`
 *   -> `{ label: 'Lum1_activity_N_P' }`
 */
export const usePlatePercentActivity = () => {
  const isCalculatingActivity = ref(false)

  const calculateActivity = async (
    plateId: number,
    settings: ActivitySettings,
  ): Promise<BackgroundCorrectionResponse> => {
    isCalculatingActivity.value = true
    try {
      return await requestApiData<BackgroundCorrectionResponse>(
        `${BACKGROUND_CORRECTION_ENDPOINT}${plateId}/${ACTIVITY_ENDPOINT}`,
        { method: 'POST', body: settings },
        ACTIVITY_ERROR_MESSAGE,
      )
    } finally {
      isCalculatingActivity.value = false
    }
  }

  return { isCalculatingActivity, calculateActivity }
}
