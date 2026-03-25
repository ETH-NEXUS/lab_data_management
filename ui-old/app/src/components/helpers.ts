import {Barcode} from 'components/models'

export type GeneralObject = {[key: string]: string}
type Segment =
  | {start: number; end: number; from: string; to: string}
  | {start: number; end: number; special: true; value: string; from: string; to: string}

const redGreenSegments: Segment[] = [
  {start: 0, end: 0.1, from: '#FF0000', to: '#ec5050'},
  {start: 0.1, end: 0.2, from: '#ec5050', to: '#CC3333'},
  {start: 0.2, end: 0.3, from: '#CC3333', to: '#993333'},
  {start: 0.3, end: 0.4, from: '#993333', to: '#660000'},
  {start: 0.4, end: 0.49, from: '#660000', to: '#330101'},
  {start: 0.49, end: 0.51, special: true, value: '#000000', from: '', to: ''},
  {start: 0.51, end: 0.6, from: '#001e00', to: '#0a510a'},
  {start: 0.6, end: 0.7, from: '#0a510a', to: '#33AA33'},
  {start: 0.7, end: 0.8, from: '#33AA33', to: '#55CC55'},
  {start: 0.8, end: 0.9, from: '#55CC55', to: '#55FF55'},
  {start: 0.9, end: 1.0, from: '#55FF55', to: '#00FF00'},
]

const greenRedSegments: Segment[] = [
  {start: 0, end: 0.1, from: '#00FF00', to: '#55FF55'},
  {start: 0.1, end: 0.2, from: '#55FF55', to: '#55CC55'},
  {start: 0.2, end: 0.3, from: '#55CC55', to: '#33AA33'},
  {start: 0.3, end: 0.4, from: '#33AA33', to: '#0a510a'},
  {start: 0.4, end: 0.49, from: '#0a510a', to: '#001e00'},
  {start: 0.49, end: 0.51, special: true, value: '#000000', from: '', to: ''},
  {start: 0.51, end: 0.6, from: '#330101', to: '#660000'},
  {start: 0.6, end: 0.7, from: '#660000', to: '#993333'},
  {start: 0.7, end: 0.8, from: '#993333', to: '#CC3333'},
  {start: 0.8, end: 0.9, from: '#CC3333', to: '#ec5050'},
  {start: 0.9, end: 1.0, from: '#ec5050', to: '#FF0000'},
]

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
  {label: 'RedGreen', value: {from: '#00FF00', to: '#FF0000'}},
]

const getSegmentColor = (p: number, segments: Segment[]): string => {
  for (const seg of segments) {
    if (p >= seg.start && p < seg.end) {
      if ('special' in seg && seg.special) {
        return seg.value
      }

      const localPct = (p - seg.start) / (seg.end - seg.start)
      return interpolateHsl(localPct, seg.from, seg.to)
    }
  }
  const last = segments[segments.length - 1]
  if ('special' in last && last.special) {
    return last.value
  }
  const localPct = (p - last.start) / (last.end - last.start)
  return interpolateHsl(localPct, last.from, last.to)
}

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
const clamp = (n: number) => Math.max(0, Math.min(1, n))

export const percentageToHsl = (
  percentage: number,
  fromColor: string,
  toColor: string,
  label: string
): string => {
  if (percentage === -1) {
    return 'rgba(255,255,255,0)'
  }

  if (label === 'RedGreen' || label === 'GreenRed') {
    const p = clamp(percentage)
    const segments = label === 'RedGreen' ? redGreenSegments : greenRedSegments
    return getSegmentColor(p, segments)
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
