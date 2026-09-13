/**
 * Threshold configuration used to mark problematic wells/plates.
 *
 * Data example:
 * - `{ id: 1, dmso: 80, amount: 2.5 }`
 */
export type Threshold = {
  id: number
  dmso: number
  amount: number
}

/**
 * API response shape returned by `GET /api/thresholds/`.
 *
 * Data example:
 * - `{ results: [{ id: 1, dmso: 80, amount: 2.5 }] }`
 */
export type ThresholdListResponse = {
  results: Threshold[]
}

/**
 * Why a well was marked: which of the two thresholds its values are below.
 */
export type ProblematicWellReason = 'volume' | 'dmso'

/**
 * One marked well with the values the instrument last reported.
 * A value is null when the instrument never reported it, which is not the
 * same as a reported zero.
 *
 * Data example:
 * - `{ position: 'I12', current_amount: 0, current_dmso: 0, reasons: ['volume', 'dmso'] }`
 */
export type ProblematicWell = {
  position: string
  current_amount: number | null
  current_dmso: number | null
  reasons: ProblematicWellReason[]
}

/**
 * Red-flag warnings grouped by library and plate barcode.
 *
 * Data example:
 * - `{
 *     "Library A": {
 *       "PLATE-001": [{ position: "I12", current_amount: 1.31, current_dmso: 94.5, reasons: ["volume"] }],
 *       "PLATE-002": []
 *     }
 *   }`
 */
export type RedFlagInfo = Record<string, Record<string, ProblematicWell[]>>

/**
 * Payload for updating threshold values.
 *
 * Data example:
 * - `{ dmso: 75, amount: 2.0 }`
 */
export type ThresholdUpdatePayload = {
  dmso: number
  amount: number
}
