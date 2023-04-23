export interface Barcode {
  NorthBarcode: string
  SouthBarcode: string
  EastBarcode: string
  WestBarcode: string
}

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

export interface ExperimentDetails {
  d: number
  project_id: number
  measurement_labels: string[]
  measurement_timestamps: {
    [label: string]: string[]
  }
  stats: {
    [label: string]: {
      [label: string]: {
        min: number[]
        max: number[]
        mean: number[]
        median: number[]
        std: number[]
        mad: number[]
      }
    }
  }
  overall_stats: {
    [label: string]: {
      min: number[]
      max: number[]
      mean: number[]
      median: number[]
      std: number[]
      mad: number[]
    }
  }
}
export interface Experiment {
  id: number
  name: string
  plates: Array<Plate>
  project: number
  available_measurement_labels: Array<string>
  details: ExperimentDetails
  barcode_specifications?: Array<BarcodeSpecification>
  description?: string
}

export interface Project {
  id: number
  name: string
  experiments: Array<Experiment>
  description?: string
  harvest_id?: number
  harvest_notes?: string
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
  templates: Array<Template>
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
  plates: Array<Plate>
}

export interface Compound {
  id: number
  name: string
  identifier: string
  structure: string
  library: CompoundLibrary
  data: object | null
  amount: number
  wells: Array<Well>
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
  hr_position?: string
  compounds?: Array<Compound>
  measurements?: Array<Measurement>
  withdrawals?: Array<Withdrawal>
  donors?: Array<Withdrawal>
  mixture?: boolean
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

export interface PlateStats {
  min: number[]
  max: number[]
  mean: number[]
  median: number[]
  std: number[]
  mad: number[]
}
export interface PlateDetails {
  id: number
  num_wells: number
  measurement_labels: Array<string>
  measurement_timestamps: {[key: string]: Array<string>}
  stats: {[key: string]: {[key: string]: PlateStats}}
  overall_stats: {[key: string]: PlateStats}
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
  compounds: Array<string>
  measurements: {[key: string]: number[]}
}

export interface Plate {
  id: number
  barcode: string
  dimension: PlateDimension
  details: PlateDetails
  wells: Array<WellDetails>
  measurement_labels?: Array<string>
  experiment?: number
  library?: number
  template?: number
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
  well: WellDetails
  position: number
}

export interface LabelValue {
  label: string
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  value: any
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

export interface ExperimentPayload {
  name?: string
  description?: string
}

export interface ProjectPayload {
  name?: string
  description?: string
}

export interface LegendColor {
  value: number
  color: string
}

export interface harvestProject {
  id: number
  value: string
  label: string
  name: string
  code: string
  is_active: boolean
  is_billable: boolean
  is_fixed_fee: boolean
  bill_by: string
  budget: null | number
  budget_by: string
  budget_is_monthly: boolean
  notify_when_over_budget: boolean
  over_budget_notification_percentage: null | number
  show_budget_to_all: boolean
  created_at: string
  updated_at: string
  starts_on: null | string
  ends_on: null | string
  over_budget_notification_date: null | string
  notes: string
  cost_budget: null | number
  cost_budget_include_expenses: null | boolean
  hourly_rate: null | number
  fee: null | number
  client: {
    id: number
    name: string
    currency: string
  }
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
}

interface DirectoryItem {
  type: 'directory'
  name: string
  children: Array<FileSystemItem>
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

export interface FormData {
  [key: string]: string | boolean
}
