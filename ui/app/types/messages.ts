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
 * Red-flag warnings grouped by library and plate barcode.
 *
 * Data example:
 * - `{
 *     "Library A": {
 *       "PLATE-001": ["A01", "B03"],
 *       "PLATE-002": ["H12"]
 *     }
 *   }`
 */
export type RedFlagInfo = Record<string, Record<string, string[]>>

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
