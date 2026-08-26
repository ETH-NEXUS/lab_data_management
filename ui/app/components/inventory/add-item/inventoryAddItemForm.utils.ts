export type AddItemFormState = {
  itemTypeId: string
  materialId: string
  roomId: string
  sectorIds: string[]
  reagentStorageTemperature: string
  reagentSafetyDataSheet: File | null
  // Optional material details. Left blank, they don't change the material at all;
  // filled in, they get PATCHed onto the selected material when the item is saved.
  additionalBrandId: string
  additionalDefaultCost: string
  additionalManufacturerId: string
  additionalVendorId: string
  additionalManufacturerCatalogNumber: string
  additionalVendorCatalogNumber: string
  additionalCapacityValue: string
  additionalCapacityUnit: string
  additionalDescription: string
  additionalSerialNumber: string
  additionalOrderNumber: string
  additionalLifetimeDays: string
  // '' = don't change, 'true' = active, 'false' = inactive.
  additionalIsActive: string
  stockUnitId: string
  quantity: string
  minimumQuantity: string
  lotNumber: string
  expiryDate: string
  notes: string
}

/**
 * Builds empty draft values for the add-item stock form.
 *
 * Returned data example:
 * - `{ itemTypeId: '', materialId: '', roomId: '', sectorIds: [], reagentStorageTemperature: '', reagentSafetyDataSheet: null, additionalBrandId: '', additionalDefaultCost: '', additionalManufacturerId: '', additionalVendorId: '', additionalManufacturerCatalogNumber: '', additionalVendorCatalogNumber: '', additionalCapacityValue: '', additionalCapacityUnit: '', additionalDescription: '', additionalSerialNumber: '', additionalOrderNumber: '', additionalLifetimeDays: '', additionalIsActive: '', stockUnitId: '', quantity: '', minimumQuantity: '', lotNumber: '', expiryDate: '', notes: '' }`
 */
export const buildInitialFormState = (): AddItemFormState => ({
  itemTypeId: '',
  materialId: '',
  roomId: '',
  sectorIds: [],
  reagentStorageTemperature: '',
  reagentSafetyDataSheet: null,
  additionalBrandId: '',
  additionalDefaultCost: '',
  additionalManufacturerId: '',
  additionalVendorId: '',
  additionalManufacturerCatalogNumber: '',
  additionalVendorCatalogNumber: '',
  additionalCapacityValue: '',
  additionalCapacityUnit: '',
  additionalDescription: '',
  additionalSerialNumber: '',
  additionalOrderNumber: '',
  additionalLifetimeDays: '',
  additionalIsActive: '',
  stockUnitId: '',
  quantity: '',
  minimumQuantity: '',
  lotNumber: '',
  expiryDate: '',
  notes: '',
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
