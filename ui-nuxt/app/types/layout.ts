import type { Plate } from '~/types/lab'

export const CONTROL_PLATES_ENDPOINT = 'plates/'
export const ADD_CONTROL_LAYOUT_ENDPOINT = 'add_control_layout/'

export const CONTROL_PLATES_PARAMS = {
  is_control_plate: 'true',
  use_as_template_to_select: 'true',
} as const

export const PROJECT_FETCH_ERROR_MESSAGE = 'Failed to load project.'
export const CONTROL_PLATES_FETCH_ERROR_MESSAGE = 'Failed to load control plates.'
export const ADD_CONTROL_LAYOUT_ERROR_MESSAGE = 'Failed to add control layout.'

export type AddControlLayoutPayload = {
  barcode_new: string
  barcode_old: string
  project_id: number | null
}

export type LayoutState = {
  controlPlates: Plate[]
}
