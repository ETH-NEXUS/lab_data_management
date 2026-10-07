import { ref } from 'vue'
import {
  BACKGROUND_CORRECTION_ENDPOINT,
  BACKGROUND_CORRECTION_ERROR_MESSAGE,
  type BackgroundCorrectionResponse,
  type BackgroundCorrectionSettings,
} from '~/types/backgroundCorrection'
import { requestApiData } from '~/utils/apiRequests'

/**
 * Starts the background correction of one plate on the server.
 *
 * Usage example:
 * - `await correctPlateBackground(42, { label: 'Lum_CTG', reference_type: 'N1', method: 'median' })`
 *   -> `{ label: 'Lum_CTG_bc_N1_median' }`
 */
export const usePlateBackgroundCorrection = () => {
  const isCorrecting = ref(false)

  const correctPlateBackground = async (
    plateId: number,
    settings: BackgroundCorrectionSettings,
  ): Promise<BackgroundCorrectionResponse> => {
    isCorrecting.value = true
    try {
      return await requestApiData<BackgroundCorrectionResponse>(
        `${BACKGROUND_CORRECTION_ENDPOINT}${plateId}/`,
        { method: 'POST', body: settings },
        BACKGROUND_CORRECTION_ERROR_MESSAGE,
      )
    } finally {
      isCorrecting.value = false
    }
  }

  return { isCorrecting, correctPlateBackground }
}
