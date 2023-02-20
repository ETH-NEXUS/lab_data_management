import {Barcode} from 'components/models'

export type GeneralObject = {[key: string]: string}

export const generateBarcodes = (prefix: string, numberOfPlates: number, sides: string[]) => {
  const shouldIncludeSide: Record<string, boolean> = {}
  for (const side of sides) {
    shouldIncludeSide[side] = true
  }

  return Array.from({length: numberOfPlates}, (_, i) => {
    const barcode = `${prefix}_${i + 1}$${prefix}_${i + 1}`
    return {
      NorthBarcode: shouldIncludeSide['North'] ? barcode : '',
      SouthBarcode: shouldIncludeSide['South'] ? barcode : '',
      EastBarcode: shouldIncludeSide['East'] ? barcode : '',
      WestBarcode: shouldIncludeSide['West'] ? barcode : '',
    }
  })
}

export const downloadCSVData = (
  columns: Array<string>,
  items: GeneralObject[] | Barcode[],
  fileName: string
): void => {
  let csv = columns.join(',')
  csv += '\r\n'

  for (const row of items) {
    csv += Object.values(row).join(',')
    csv += '\r\n'
  }

  const blob = new Blob([csv], {type: 'text/csv;charset=utf-8;'})
  const link = document.createElement('a')
  link.href = URL.createObjectURL(blob)
  link.download = fileName
  link.click()
}
