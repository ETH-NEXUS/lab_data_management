import type { BarcodeSide } from '~/types/experiments'
import type { Barcode } from '~/types/lab'

/**
 * Select options for barcode side selection in forms.
 *
 * Data example:
 * - `[{ label: 'North', value: 'North' }, { label: 'South', value: 'South' }]`
 */
export const BARCODE_SIDE_OPTIONS: Array<{ label: BarcodeSide; value: BarcodeSide }> = [
  { label: 'North', value: 'North' },
  { label: 'South', value: 'South' },
  { label: 'East', value: 'East' },
  { label: 'West', value: 'West' },
]

/**
 * Stable CSV column order for generated barcode exports.
 *
 * Data example:
 * - `['NorthBarcode', 'SouthBarcode', 'EastBarcode', 'WestBarcode']`
 */
export const BARCODE_CSV_COLUMNS: Array<keyof Barcode> = ['NorthBarcode', 'SouthBarcode', 'EastBarcode', 'WestBarcode']

/**
 * Generates barcode rows from one specification.
 *
 * Accepted input example:
 * - `{ prefix: 'A001', numberOfPlates: 2, sides: ['North', 'West'] }`
 *
 * Returned data example:
 * - `[
 *      { NorthBarcode: 'A001_1$A001_1', SouthBarcode: '', EastBarcode: '', WestBarcode: 'A001_1$A001_1' },
 *      { NorthBarcode: 'A001_2$A001_2', SouthBarcode: '', EastBarcode: '', WestBarcode: 'A001_2$A001_2' }
 *    ]`
 */
export const generateBarcodes = (prefix: string, numberOfPlates: number, sides: string[]): Barcode[] => {
  const includedSides: Record<string, boolean> = {}
  for (const side of sides) {
    includedSides[side] = true
  }

  return Array.from({ length: numberOfPlates }, (_, index) => {
    const barcodeValue = `${prefix}_${index + 1}$${prefix}_${index + 1}`
    return {
      NorthBarcode: includedSides.North ? barcodeValue : '',
      SouthBarcode: includedSides.South ? barcodeValue : '',
      EastBarcode: includedSides.East ? barcodeValue : '',
      WestBarcode: includedSides.West ? barcodeValue : '',
    }
  })
}

/**
 * Escapes one CSV field value.
 *
 * Accepted value example:
 * - `North,plate`
 *
 * Returned value example:
 * - `"North,plate"`
 */
const escapeCsvValue = (value: unknown): string => {
  const text = String(value ?? '')
  const escapedText = text.replace(/"/g, '""')
  const needsQuotes = /[",\n\r]/.test(escapedText)
  if (!needsQuotes) return escapedText
  return `"${escapedText}"`
}

/**
 * Downloads CSV data in the browser.
 *
 * Accepted data example:
 * - `columns = ['NorthBarcode', 'SouthBarcode']`
 * - `rows = [{ NorthBarcode: 'A001_1$A001_1', SouthBarcode: '' }]`
 * - `fileName = 'barcodes.csv'`
 */
export const downloadCsvData = <TRow extends Record<string, unknown>>(
  columns: string[],
  rows: TRow[],
  fileName: string,
): void => {
  const headerLine = columns.map((column) => escapeCsvValue(column)).join(',')
  const dataLines = rows.map((row) => columns.map((column) => escapeCsvValue(row[column])).join(','))
  const csvContent = [headerLine, ...dataLines].join('\r\n')

  const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' })
  const url = URL.createObjectURL(blob)
  const link = document.createElement('a')
  link.href = url
  link.download = fileName
  document.body.appendChild(link)
  link.click()
  document.body.removeChild(link)
  URL.revokeObjectURL(url)
}
