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

export type Palette = {
  label: string
  value: {
    from: string
    to: string
  }
}

export const palettes: Array<Palette> = [
  {label: 'OrangeRed', value: {from: '#fff7bc', to: '#993404'}},
  {label: 'GreenRed', value: {from: '#FF0000', to: '#00FF00'}},
  {label: 'Green', value: {from: '#c7e9c0', to: '#006d2c'}},
  {label: 'Blue', value: {from: '#92c5de', to: '#0b2746'}},
  {label: 'GreenBrown', value: {from: '#b8e186', to: '#662506'}},
  {label: 'Grey', value: {from: '#ffffff', to: '#525252'}},
  {label: 'Magma', value: {from: '#ffffff', to: '#000000'}},
]

function easeInOutCubic(t: number) {
  return t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2
}

function interpolateHsl(p: number, colorA: string, colorB: string) {
  const from = hexToRgb(colorA)
  const to = hexToRgb(colorB)
  const easedT = easeInOutCubic(p)

  const hue = easedT * (to.h - from.h) + from.h
  const sat = easedT * (to.s - from.s) + from.s
  const lig = easedT * (to.l - from.l) + from.l
  return `hsl(${hue}, ${sat}%, ${lig}%)`
}

export const percentageToHsl = (percentage: number, fromColor: string, toColor: string) => {
  // if percentage is not given (-1) we return a transparent color
  if (percentage === -1) {
    return 'rgba(255,255,255,0)'
  }

  // SPECIAL CASE: GreenRed palettewe  want black in the middle

  if (fromColor.toLowerCase() === '#ff0000' && toColor.toLowerCase() === '#00ff00') {
    const p = Math.max(0, Math.min(1, percentage)) // p в пределах [0..1]

    // ----- 0..10% -----
    if (p < 0.1) {
      const localPct = p / 0.1
      return interpolateHsl(localPct, '#00FF00', '#55FF55')
    }
    // ----- 10..20% -----
    else if (p < 0.2) {
      const localPct = (p - 0.1) / 0.1
      return interpolateHsl(localPct, '#55FF55', '#55CC55')
    }
    // ----- 20..30% -----
    else if (p < 0.3) {
      const localPct = (p - 0.2) / 0.1
      return interpolateHsl(localPct, '#55CC55', '#33AA33')
    }
    // ----- 30..40% -----
    else if (p < 0.4) {
      const localPct = (p - 0.3) / 0.1
      return interpolateHsl(localPct, '#33AA33', '#147514')
    }
    // ----- 40..49% (9%) -----
    else if (p < 0.49) {
      const localPct = (p - 0.4) / 0.09
      return interpolateHsl(localPct, '#147514', '#001e00')
    }
    // ----- 49..51% →
    else if (p < 0.51) {
      return '#000000'
    }
    // ----- 51..60% (9%) -----
    else if (p < 0.6) {
      const localPct = (p - 0.51) / 0.09
      return interpolateHsl(localPct, '#330101', '#660000')
    }
    // ----- 60..70% -----
    else if (p < 0.7) {
      const localPct = (p - 0.6) / 0.1
      return interpolateHsl(localPct, '#660000', '#993333')
    }
    // ----- 70..80% -----
    else if (p < 0.8) {
      const localPct = (p - 0.7) / 0.1
      return interpolateHsl(localPct, '#993333', '#CC3333')
    }
    // ----- 80..90% -----
    else if (p < 0.9) {
      const localPct = (p - 0.8) / 0.1
      return interpolateHsl(localPct, '#CC3333', '#ec5050')
    }
    // ----- 90..100% -----
    else {
      const localPct = (p - 0.9) / 0.1
      return interpolateHsl(localPct, '#ec5050', '#FF0000')
    }
  }

  //  existing single-step interpolation for all other palettes
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
