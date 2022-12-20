import {CompoundExpressionNode} from '@vue/compiler-core'

export interface Todo {
    id: number
    content: string
}

export interface Project {
    name: string
}

export interface Meta {
    totalCount: number
}

export interface Sample {
    name: string
}

export interface CompoundLibrary {
    name: string
    file_name: string
    plates: Array<Plate>
}

export interface Compound {
    name: string
    identifier: string
    structure: string
    library: CompoundLibrary
    data: object | null
}

export interface Well {
    id: number
    plate: Plate
    position: number
    hr_position: string
    sample: Sample
    compounds: Array<Compound>
    amount: number
    measurements: Array<Measurement>
}

export interface PlateDimension {
    name: string
    cols: number
    rows: number
}

export interface Plate {
    barcode: string
    dimension: PlateDimension
    project?: Project
    library?: CompoundLibrary
    wells?: Array<Well>
}

export interface Measurement {
    name: string
    abbrev: string
    unit: string
    value: number
}
