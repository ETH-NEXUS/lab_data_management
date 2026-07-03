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
const parseInteger = (value: string | number | null | undefined): number | null => {
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
