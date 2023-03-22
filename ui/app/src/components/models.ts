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

export interface Experiment {
  id: number
  name: string
  plates: Array<Plate>
  project: number
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
  hr_position?: string
  compounds?: Array<Compound>
  measurements?: Array<Measurement>
  withdrawals?: Array<Withdrawal>
  donors?: Array<Withdrawal>
  mixture?: boolean
  status: string
}

export interface PlateDimension {
  id: number
  name: string
  cols: number
  rows: number
}

export interface MinMax {
  min: number
  max: number
  min_all_types: number
  max_all_types: number
}

export interface Plate {
  id: number
  barcode: string
  dimension: PlateDimension
  measurements: Array<string>
  z_primes: {[key: string]: number}
  min_max: {[key: string]: MinMax}
  experiment?: Experiment
  library?: CompoundLibrary
  template?: Template
  wells?: Array<Well>
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
  feature: MeasurementFeature
  well?: Well
}

export interface WellInfo {
  well: Well
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
