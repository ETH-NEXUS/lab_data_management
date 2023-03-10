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

export const percentageToHsl = (percentage: number, fromColor: string, toColor: string) => {
  // if percentage is not given (-1) we return a transparent color
  if (percentage === -1) {
    return 'rgba(255,255,255,0)'
  }

  const fromRgb = hexToRgb(fromColor)
  const toRgb = hexToRgb(toColor)

  const hue = percentage * (toRgb.h - fromRgb.h) + fromRgb.h
  const saturation = percentage * (toRgb.s - fromRgb.s) + fromRgb.s
  const lightness = percentage * (toRgb.l - fromRgb.l) + fromRgb.l

  return `hsl(${hue}, ${saturation}%, ${lightness}%)`
}

export const hexToRgb = (hex: string) => {
  const r = parseInt(hex.slice(1, 3), 16)
  const g = parseInt(hex.slice(3, 5), 16)
  const b = parseInt(hex.slice(5, 7), 16)
  const {h, s, l} = rgbToHsl(r, g, b)
  return {r, g, b, h, s, l}
}

const rgbToHsl = (r: number, g: number, b: number) => {
  r /= 255
  g /= 255
  b /= 255

  const max = Math.max(r, g, b),
    min = Math.min(r, g, b)
  let h = (max + min) / 2
  let s = (max + min) / 2
  const l = (max + min) / 2

  if (max == min) {
    h = s = 0
  } else {
    const d = max - min
    s = l > 0.5 ? d / (2 - max - min) : d / (max + min)
    switch (max) {
      case r:
        h = (g - b) / d + (g < b ? 6 : 0)
        break
      case g:
        h = (b - r) / d + 2
        break
      case b:
        h = (r - g) / d + 4
        break
    }
    if (h) {
      h /= 6
    }
  }

  return {h: h * 360, s: s * 100, l: l * 100}
}
