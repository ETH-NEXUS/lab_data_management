export type AddItemFormState = {
  itemTypeId: string
  materialId: string
  roomId: string
  sectorId: string
  stockUnitId: string
  quantity: string
  minimumQuantity: string
  lotNumber: string
  expiryDate: string
  notes: string
  isFavorite: boolean
}

/**
 * Builds empty draft values for the add-item stock form.
 *
 * Returned data example:
 * - `{ itemTypeId: '', materialId: '', roomId: '', sectorId: '', stockUnitId: '', quantity: '', minimumQuantity: '', lotNumber: '', expiryDate: '', notes: '', isFavorite: false }`
 */
export const buildInitialFormState = (): AddItemFormState => ({
  itemTypeId: '',
  materialId: '',
  roomId: '',
  sectorId: '',
  stockUnitId: '',
  quantity: '',
  minimumQuantity: '',
  lotNumber: '',
  expiryDate: '',
  notes: '',
  isFavorite: false,
})

/**
 * Parses decimal-like text into a finite number.
 *
 * Input examples:
 * - `'2.5'`
 * - `'96'`
 * - `null`
 *
 * Returned examples:
 * - `2.5`
 * - `96`
 * - `null`
 */
export const parseDecimal = (value: string | null | undefined): number | null => {
  if (value == null) {
    return null
  }

  const parsedValue = Number.parseFloat(value.trim())
  if (!Number.isFinite(parsedValue)) {
    return null
  }

  return parsedValue
}

/**
 * Formats decimal values for API-compatible string fields.
 *
 * Input examples:
 * - `2`
 * - `2.041666666`
 *
 * Returned examples:
 * - `'2'`
 * - `'2.041667'`
 */
export const formatDecimal = (value: number): string => {
  const fixedValue = value.toFixed(6)
  const normalizedValue = fixedValue.replace(/\.?0+$/, '')
  return normalizedValue === '' ? '0' : normalizedValue
}
