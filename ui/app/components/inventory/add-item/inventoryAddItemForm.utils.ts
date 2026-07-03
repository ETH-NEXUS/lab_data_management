export type AddItemFormState = {
  itemTypeId: string
  materialId: string
  roomId: string
  sectorIds: string[]
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
 * - `{ itemTypeId: '', materialId: '', roomId: '', sectorIds: [], stockUnitId: '', quantity: '', minimumQuantity: '', lotNumber: '', expiryDate: '', notes: '', isFavorite: false }`
 */
export const buildInitialFormState = (): AddItemFormState => ({
  itemTypeId: '',
  materialId: '',
  roomId: '',
  sectorIds: [],
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
 * Parses integer-like text into a finite whole number.
 *
 * Input examples:
 * - `'55'`
 * - `55`
 * - `'0'`
 * - `'2.5'`
 * - `null`
 *
 * Returned examples:
 * - `55`
 * - `55`
 * - `0`
 * - `null`
 * - `null`
 */
export const parseInteger = (value: string | number | null | undefined): number | null => {
  if (value == null) {
    return null
  }

  if (typeof value === 'number') {
    return Number.isInteger(value) && value >= 0 ? value : null
  }

  const normalizedValue = value.trim()
  if (normalizedValue === '') {
    return null
  }

  const parsedValue = Number(normalizedValue)
  if (!Number.isInteger(parsedValue) || parsedValue < 0) {
    return null
  }

  return parsedValue
}

/**
 * Parses a list of integer-like values into positive integer identifiers.
 *
 * Input examples:
 * - `['5', '9']`
 * - `['5', 'bad', '9']`
 * - `[]`
 *
 * Returned examples:
 * - `[5, 9]`
 * - `[5, 9]`
 * - `[]`
 */
export const parsePositiveIntegerList = (values: string[] | null | undefined): number[] => {
  if (!Array.isArray(values)) {
    return []
  }

  return values.reduce<number[]>((parsedValues, value) => {
    const parsedValue = parseInteger(value)
    if (parsedValue === null || parsedValue <= 0) {
      return parsedValues
    }

    parsedValues.push(parsedValue)
    return parsedValues
  }, [])
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
