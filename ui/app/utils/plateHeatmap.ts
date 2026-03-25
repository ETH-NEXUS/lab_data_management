import type { PlatePaletteOption } from '~/types/plates'
import type { LegendColor } from '~/types/lab'

/**
 * Default heatmap palette options migrated from the old UI helper.
 *
 * Data examples:
 * - `{ label: 'GreenRed', value: { from: '#FF0000', to: '#00FF00' } }`
 * - `{ label: 'Blue', value: { from: '#92c5de', to: '#0b2746' } }`
 */
export const platePalettes: PlatePaletteOption[] = [
  { label: 'OrangeRed', value: { from: '#fff7bc', to: '#993404' } },
  { label: 'GreenRed', value: { from: '#FF0000', to: '#00FF00' } },
  { label: 'Green', value: { from: '#c7e9c0', to: '#006d2c' } },
  { label: 'Blue', value: { from: '#92c5de', to: '#0b2746' } },
  { label: 'GreenBrown', value: { from: '#b8e186', to: '#662506' } },
  { label: 'Grey', value: { from: '#ffffff', to: '#525252' } },
  { label: 'RedGreen', value: { from: '#00FF00', to: '#FF0000' } },
]

/**
 * Returns the preferred default palette.
 *
 * Returned data example:
 * - `{ label: 'GreenRed', value: { from: '#FF0000', to: '#00FF00' } }`
 */
export const getDefaultPlatePalette = (): PlatePaletteOption => {
  return platePalettes.find((palette) => palette.label === 'GreenRed') ?? platePalettes[0]!
}

const clamp = (value: number, min: number, max: number): number => Math.min(max, Math.max(min, value))

type Segment =
  | { start: number; end: number; from: string; to: string }
  | { start: number; end: number; special: true; value: string; from: string; to: string }

const redGreenSegments: Segment[] = [
  { start: 0, end: 0.1, from: '#FF0000', to: '#ec5050' },
  { start: 0.1, end: 0.2, from: '#ec5050', to: '#CC3333' },
  { start: 0.2, end: 0.3, from: '#CC3333', to: '#993333' },
  { start: 0.3, end: 0.4, from: '#993333', to: '#660000' },
  { start: 0.4, end: 0.49, from: '#660000', to: '#330101' },
  { start: 0.49, end: 0.51, special: true, value: '#000000', from: '', to: '' },
  { start: 0.51, end: 0.6, from: '#001e00', to: '#0a510a' },
  { start: 0.6, end: 0.7, from: '#0a510a', to: '#33AA33' },
  { start: 0.7, end: 0.8, from: '#33AA33', to: '#55CC55' },
  { start: 0.8, end: 0.9, from: '#55CC55', to: '#55FF55' },
  { start: 0.9, end: 1.0, from: '#55FF55', to: '#00FF00' },
]

const greenRedSegments: Segment[] = [
  { start: 0, end: 0.1, from: '#00FF00', to: '#55FF55' },
  { start: 0.1, end: 0.2, from: '#55FF55', to: '#55CC55' },
  { start: 0.2, end: 0.3, from: '#55CC55', to: '#33AA33' },
  { start: 0.3, end: 0.4, from: '#33AA33', to: '#0a510a' },
  { start: 0.4, end: 0.49, from: '#0a510a', to: '#001e00' },
  { start: 0.49, end: 0.51, special: true, value: '#000000', from: '', to: '' },
  { start: 0.51, end: 0.6, from: '#330101', to: '#660000' },
  { start: 0.6, end: 0.7, from: '#660000', to: '#993333' },
  { start: 0.7, end: 0.8, from: '#993333', to: '#CC3333' },
  { start: 0.8, end: 0.9, from: '#CC3333', to: '#ec5050' },
  { start: 0.9, end: 1.0, from: '#ec5050', to: '#FF0000' },
]

const hexToRgb = (hex: string): { r: number; g: number; b: number } => {
  const normalizedHex = hex.replace('#', '')

  if (normalizedHex.length !== 6) {
    return { r: 0, g: 0, b: 0 }
  }

  return {
    r: Number.parseInt(normalizedHex.slice(0, 2), 16),
    g: Number.parseInt(normalizedHex.slice(2, 4), 16),
    b: Number.parseInt(normalizedHex.slice(4, 6), 16),
  }
}

const rgbToHex = (r: number, g: number, b: number): string => {
  const toHex = (value: number): string => value.toString(16).padStart(2, '0')
  return `#${toHex(clamp(Math.round(r), 0, 255))}${toHex(clamp(Math.round(g), 0, 255))}${toHex(clamp(Math.round(b), 0, 255))}`
}

const rgbToHsl = (r: number, g: number, b: number) => {
  r /= 255
  g /= 255
  b /= 255

  const max = Math.max(r, g, b)
  const min = Math.min(r, g, b)
  let h = (max + min) / 2
  let s = (max + min) / 2
  const l = (max + min) / 2

  if (max === min) {
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
    h /= 6
  }

  return { h: h * 360, s: s * 100, l: l * 100 }
}

const hexToHsl = (hex: string) => {
  const rgb = hexToRgb(hex)
  return rgbToHsl(rgb.r, rgb.g, rgb.b)
}

const easeInOutCubic = (t: number) => (t < 0.5 ? 4 * t * t * t : 1 - (-2 * t + 2) ** 3 / 2)

const interpolateHsl = (p: number, colorA: string, colorB: string) => {
  const from = hexToHsl(colorA)
  const to = hexToHsl(colorB)
  const easedT = easeInOutCubic(p)

  const hue = easedT * (to.h - from.h) + from.h
  const sat = easedT * (to.s - from.s) + from.s
  const lig = easedT * (to.l - from.l) + from.l
  return `hsl(${hue}, ${sat}%, ${lig}%)`
}

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
  if (!last) return 'transparent'

  if ('special' in last && last.special) {
    return last.value
  }
  const localPct = (p - last.start) / (last.end - last.start)
  return interpolateHsl(localPct, last.from, last.to)
}

/**
 * Interpolates between two palette colors with linear percentage.
 *
 * Accepted data examples:
 * - `from = '#FF0000', to = '#00FF00', percent = 0.5`
 *
 * Returned data examples:
 * - `'#808000'`
 */
export const interpolatePaletteColor = (from: string, to: string, percent: number): string => {
  const safePercent = clamp(percent, 0, 1)
  const fromRgb = hexToRgb(from)
  const toRgb = hexToRgb(to)

  const r = fromRgb.r + (toRgb.r - fromRgb.r) * safePercent
  const g = fromRgb.g + (toRgb.g - fromRgb.g) * safePercent
  const b = fromRgb.b + (toRgb.b - fromRgb.b) * safePercent

  return rgbToHex(r, g, b)
}

/**
 * Legacy-compatible heatmap color interpolation by percentage.
 *
 * Accepted data examples:
 * - `percentage = 0.55, fromColor = '#FF0000', toColor = '#00FF00', label = 'GreenRed'`
 *
 * Returned data examples:
 * - `'hsl(120, 100%, 50%)'`
 * - `'rgba(255,255,255,0)'` when percentage is invalid marker `-1`
 */
export const percentageToHsl = (percentageValue: number, fromColor: string, toColor: string, label: string): string => {
  if (percentageValue === -1) {
    return 'rgba(255,255,255,0)'
  }

  if (label === 'RedGreen' || label === 'GreenRed') {
    const p = clamp(percentageValue, 0, 1)
    const segments = label === 'RedGreen' ? redGreenSegments : greenRedSegments
    return getSegmentColor(p, segments)
  }

  const fromHsl = hexToHsl(fromColor)
  const toHsl = hexToHsl(toColor)
  const hue = percentageValue * (toHsl.h - fromHsl.h) + fromHsl.h
  const saturation = percentageValue * (toHsl.s - fromHsl.s) + fromHsl.s
  const lightness = percentageValue * (toHsl.l - fromHsl.l) + fromHsl.l
  return `hsl(${hue}, ${saturation}%, ${lightness}%)`
}

/**
 * Converts one measurement value into a heatmap color.
 *
 * Accepted data examples:
 * - `value = 220, min = 0, max = 1000, palette = { label: 'GreenRed', value: { from: '#FF0000', to: '#00FF00' } }`
 *
 * Returned data examples:
 * - `'#c73800'`
 * - `'transparent'` for invalid ranges
 */
export const getHeatmapColor = (
  value: number | null,
  min: number,
  max: number,
  palette: PlatePaletteOption,
): string => {
  if (value === null || Number.isNaN(value) || Number.isNaN(min) || Number.isNaN(max) || min === max) {
    return 'transparent'
  }

  const percent = (value - min) / (max - min)
  return percentageToHsl(percent, palette.value.from, palette.value.to, palette.label)
}

/**
 * Creates the stepped plate legend used by the old UI.
 *
 * Returned data example:
 * - `[{ value: 1.2, color: 'hsl(...)' }, ...]` (reversed order for top-down display)
 */
export const buildPlateLegend = (
  min: number,
  max: number,
  palette: PlatePaletteOption,
  numberOfSteps = 20,
): LegendColor[] => {
  const legend: LegendColor[] = []
  const step = (max - min) / numberOfSteps

  for (let i = 0; i <= numberOfSteps; i++) {
    const value = min + i * step
    const color = percentageToHsl((value - min) / (max - min), palette.value.from, palette.value.to, palette.label)
    legend.push({ value, color })
  }

  return legend.reverse()
}
