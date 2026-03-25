/**
 * Shared domain types migrated from the legacy UI.
 *
 * These are intentionally grouped in one thematic file because many of them
 * reference each other (for example `Well`, `Compound`, `Withdrawal`, `Plate`).
 * Keeping them together avoids circular-import issues and keeps usage simple.
 */

/**
 * Barcode generated for a plate side combination.
 *
 * Data example:
 * - `{ NorthBarcode: 'A001-N', SouthBarcode: 'A001-S', EastBarcode: 'A001-E', WestBarcode: 'A001-W' }`
 */
export interface Barcode {
  NorthBarcode: string
  SouthBarcode: string
  EastBarcode: string
  WestBarcode: string
}

/**
 * Barcode specification configuration attached to an experiment.
 *
 * Data example:
 * - `{ id: 1, prefix: 'A001', number_of_plates: 12, sides: ['North', 'South'], experiment: 21 }`
 */
export interface BarcodeSpecification {
  id: number
  created_at: Date | string
  modified_at: Date | string
  prefix: string
  number_of_plates: number
  sides: string[]
  experiment: number
}

export interface Todo {
  id: number
  content: string
}

/**
 * Stats shape used by plate/experiment details.
 *
 * Data example:
 * - `{ min: [0.1], max: [0.9], mean: [0.4], median: [0.4], std: [0.05], mad: [0.03] }`
 */
export interface PlateStats {
  min: number[]
  max: number[]
  mean: number[]
  median: number[]
  std: number[]
  mad: number[]
}

/**
 * Experiment-level computed details.
 *
 * Data example:
 * - `{ d: 384, project_id: 5, measurement_labels: ['OD600'], measurement_timestamps: { OD600: ['2026-01-01T00:00:00Z'] }, stats: {}, overall_stats: {} }`
 */
export interface ExperimentDetails {
  d: number
  project_id: number
  measurement_labels: string[]
  measurement_timestamps: {
    [label: string]: string[]
  }
  stats: {
    [label: string]: {
      [label: string]: PlateStats
    }
  }
  overall_stats: {
    [label: string]: PlateStats
  }
}

export interface Template {
  id: number
  name: string
  plate: Plate
  category: number
}

export interface TemplateCategory {
  id: number
  name: string
  templates: Template[]
}

export interface Meta {
  totalCount: number
}

export interface Sample {
  id: number
  name: string
}

export interface CompoundLibrary {
  id: number
  name: string
  file_name: string
  plates: Plate[]
}

export interface Compound {
  id: number
  name: string
  identifier: string
  structure: string
  library: CompoundLibrary
  data: object | null
  amount: number
  wells: Well[]
  compound: number
}

export interface Withdrawal {
  id: number
  created_at: Date
  well: Well
  amount: number
  target_well?: Well
}

export interface Well {
  id: number
  plate: Plate
  position: number
  sample: Sample
  amount: number
  type: string
  status: string
  is_invalid: boolean
  hr_position?: string
  compounds?: Compound[]
  measurements?: Measurement[]
  withdrawals?: Withdrawal[]
  donors?: Withdrawal[]
  mixture?: boolean
  current_info?: {
    current_amount: number
    current_dmso: number
  }
}

export interface PlateDimension {
  id: number
  name: string
  cols: number
  rows: number
}

export interface MinMax {
  feature: string
  timestamp: string
  min: number
  max: number
  min_all_types: number
  max_all_types: number
}

export interface PlateDetails {
  id: number
  num_wells: number
  measurement_labels: string[]
  measurement_timestamps: { [key: string]: string[] }
  stats: { [key: string]: { [key: string]: PlateStats } }
  overall_stats: { [key: string]: PlateStats }
}

export interface WellDetails {
  id: number
  plate_id: number
  type: string
  status: string
  position: number
  hr_position: string
  initial_amount: number
  withdrawal: number
  amount: number
  compounds: string[]
  measurements: { [key: string]: number[] }
}

export interface Plate {
  id: number
  barcode: string
  dimension: PlateDimension
  details: PlateDetails
  wells: WellDetails[]
  measurement_labels?: string[]
  experiment?: number | IdName | null
  library?: number | IdName | null
  template?: number | IdName | null
  is_control_plate?: boolean
  use_as_template_to_select?: boolean
  archived?: boolean
  status?: string
}

export interface MeasurementFeature {
  id: number
  name: string
  abbrev: string
  unit: string
}

export interface Measurement {
  id: number
  value: number
  label: string
  feature: MeasurementFeature
  measured_at: Date | string
  well?: Well
}

export interface WellInfo {
  well: WellDetails | undefined
  position: number
}

export interface LabelValue {
  label: string
  value: unknown
}

export interface IdName {
  id: number
  name: string
}

export interface PlateLabelValue {
  label: string
  value: number
  library: IdName | null
  experiment: IdName | null
}

export interface PlateMapping {
  source_plate: number | undefined
  target_plate: number | undefined
  from_column: string | undefined
  to_column: string | undefined
  amount_column: string | undefined
  delimiter: string | undefined
  quotechar: string | undefined
  mapping_file: File | undefined
  amount: number | undefined
}

export interface DimensionsOption {
  label: string
  value: number
}

export interface LegendColor {
  value: number
  color: string
}

export interface SelectOption<T> {
  label: string
  value: T
}

export interface TimeSeriesChart {
  series: {
    name?: string
    data?: string[] | number[]
  }[]
  colors: string[]
  chart: {
    id: string
    type: string
    height: number
    zoom: {
      enabled: boolean
    }
    dropShadow: {
      enabled: boolean
      color: string
      top: number
      left: number
      blur: number
      opacity: number
    }
  }
  dataLabels: {
    enabled: boolean
  }
  stroke: {
    curve: string
  }
  title: {
    text?: string
    align: string
  }
  grid: {
    borderColor: string
    row: {
      colors: string[]
      opacity: number
    }
  }
  xaxis: {
    categories?: string[] | Date[] | number[] | (string | Date | number)[] | undefined
    type: string
  }
  yaxis: {
    title: {
      text: string
    }
  }
  markers: {
    size: number
  }
  legend: {
    position: string
    horizontalAlign: string
    floating: boolean
    offsetY: number
    offsetX: number
  }
}

export interface PlotData {
  name: string
  data: string[] | number[] | Date[] | (string | Date | number)[]
  color: string
  categories: string[] | number[] | Date[] | (string | Date | number)[]
}

export interface CalculatorPayload {
  expression?: string
  new_label?: string
  used_labels?: string[]
  separate_time_series_points?: boolean
  plate_id?: number | null
  experiment_id?: number | null
}

interface FileItem {
  type: 'file'
  name: string
  path: string
}

interface DirectoryItem {
  type: 'directory'
  name: string
  path: string
  children: FileSystemItem[]
}

export type FileSystemItem = FileItem | DirectoryItem

export interface Options {
  [key: string]: {
    type: 'str' | 'bool'
    label: string
    required: boolean
    inputType?: string
    choices?: string[]
  }
}

export interface GeneralFormData {
  [key: string]: string | boolean
}

/**
 * Plate-level metadata used by report/prefill workflows.
 *
 * Data example:
 * - `{ measurement_label: 'OD600', measurement_timestamp: ['2026-01-01T00:00:00Z'], plate_barcode: 'A001', lib_plate_barcode: 'LIB-1', replicate: '1', cell_type: 'HeLa', condition: 'Control' }`
 */
export interface PlateInfo {
  measurement_label: string
  measurement_timestamp: string[]
  plate_barcode: string
  lib_plate_barcode: string
  replicate: string
  cell_type: string
  condition: string
}
