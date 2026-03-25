import type { Options } from '~/types/lab'

export type ManagementTabKey = 'map' | 'importControlPlate' | 'importLibraryPlate' | 'importSdf'

export type ManagementCommandConfig = {
  tab: ManagementTabKey
  command: string
  what: string
  labelKey: string
  options: Options
}

/**
 * Form fields for the "map" command.
 *
 * Accepted data example:
 * - `{ machine: 'echo', path: '/data/mappings', mapping_file: 'headers.yml', experiment_name: 'Run 2026-03', measurement_name: 'Label1' }`
 */
export const mapCommandOptions: Options = {
  machine: {
    type: 'str',
    choices: ['echo', 'm1000', 'C10-imager', 'C10-reader'],
    label: 'Machine to map from',
    required: true,
  },
  path: {
    type: 'str',
    label: 'Path to the directory containing the mapping files',
    required: true,
  },
  mapping_file: {
    type: 'str',
    label: 'For echo mapping only: YML file with column headers. If empty, default headers are used.',
    required: false,
  },
  experiment_name: {
    type: 'str',
    label: 'Experiment name',
    required: true,
  },
  measurement_name: {
    type: 'str',
    label: 'Measurement name. If empty, default "Label1" is used.',
    required: false,
  },
}

/**
 * Form fields for "import control plate".
 *
 * Accepted data example:
 * - `{ input_file: '/data/imports/control.csv', project_name: 'Project A', plate_barcode: 'CTRL-001' }`
 */
export const importControlPlateOptions: Options = {
  input_file: {
    type: 'str',
    label: 'Input file path',
    required: true,
  },
  project_name: {
    type: 'str',
    label: 'Project name',
    required: true,
  },
  plate_barcode: {
    type: 'str',
    label: 'Plate barcode',
    required: false,
  },
}

/**
 * Form fields for "import library plate".
 *
 * Accepted data example:
 * - `{ input_file: '/data/imports/library.csv', library_name: 'Kinase', plate_barcode: 'LIB-001' }`
 */
export const importLibraryPlateOptions: Options = {
  input_file: {
    type: 'str',
    label: 'Input file path',
    required: true,
  },
  library_name: {
    type: 'str',
    label: 'Library name',
    required: false,
  },
  plate_barcode: {
    type: 'str',
    label: 'Plate barcode',
    required: false,
  },
}

/**
 * Form fields for "import sdf".
 *
 * Accepted data example:
 * - `{ input_file: '/data/imports/library.sdf', library_name: 'Chem Collection', mapping_file: '/data/imports/sdf_mapping.yml' }`
 */
export const importSdfOptions: Options = {
  input_file: {
    type: 'str',
    label: 'Input file path',
    required: true,
  },
  library_name: {
    type: 'str',
    label: 'Library name',
    required: false,
  },
  mapping_file: {
    type: 'str',
    label: 'Mapping file for SDF columns. If empty, default mapping is used.',
    required: false,
  },
}

/**
 * Tab definitions used by the migrated management page.
 *
 * Data example:
 * - `[{ tab: 'map', command: 'map', what: '', labelKey: 'management.echo' }]`
 */
export const managementCommandConfigs: ManagementCommandConfig[] = [
  {
    tab: 'map',
    command: 'map',
    what: '',
    labelKey: 'management.echo',
    options: mapCommandOptions,
  },
  {
    tab: 'importControlPlate',
    command: 'import',
    what: 'control_plate',
    labelKey: 'management.import_control_plate',
    options: importControlPlateOptions,
  },
  {
    tab: 'importLibraryPlate',
    command: 'import',
    what: 'library_plate',
    labelKey: 'management.import_library_plate',
    options: importLibraryPlateOptions,
  },
  {
    tab: 'importSdf',
    command: 'import',
    what: 'sdf',
    labelKey: 'management.import_sdf_library',
    options: importSdfOptions,
  },
]
