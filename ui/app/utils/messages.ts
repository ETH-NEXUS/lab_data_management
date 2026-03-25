import type { Threshold } from '~/types/messages'

/**
 * API endpoints used by the messages page.
 *
 * Values example:
 * - `THRESHOLDS_ENDPOINT = 'thresholds/'`
 * - `RED_FLAG_ENDPOINT = 'compoundlib/redflag/'`
 */
export const THRESHOLDS_ENDPOINT = 'thresholds/'
export const RED_FLAG_ENDPOINT = 'compoundlib/redflag/'
export const RECALCULATE_STATUS_ENDPOINT = 'compoundlib/recalculate_status/'

/**
 * Fallback threshold used before API data is loaded.
 *
 * Data example:
 * - `{ id: 1, dmso: 80, amount: 2.5 }`
 */
export const DEFAULT_THRESHOLD: Threshold = {
  id: 1,
  dmso: 80,
  amount: 2.5,
}

/**
 * Normalizes unknown threshold-like input into a safe `Threshold` object.
 *
 * Accepted input examples:
 * - `null`
 * - `{ id: 2, dmso: '70', amount: '1.5' }`
 *
 * Returned data example:
 * - `{ id: 2, dmso: 70, amount: 1.5 }`
 */
export const normalizeThreshold = (value: Partial<Threshold> | null | undefined): Threshold => {
  const id = Number(value?.id ?? DEFAULT_THRESHOLD.id)
  const dmso = Number(value?.dmso ?? DEFAULT_THRESHOLD.dmso)
  const amount = Number(value?.amount ?? DEFAULT_THRESHOLD.amount)

  return {
    id: Number.isNaN(id) ? DEFAULT_THRESHOLD.id : id,
    dmso: Number.isNaN(dmso) ? DEFAULT_THRESHOLD.dmso : dmso,
    amount: Number.isNaN(amount) ? DEFAULT_THRESHOLD.amount : amount,
  }
}
