import {Options} from 'components/models'

interface SideData {
  id: string
  label: string
  styleCheckbox: string
  styleLabel: string
}

export const sidesData: SideData[] = [
  {
    id: 'top-left',
    label: 'North',
    styleCheckbox: 'position: absolute; top: -20px; left: 35px',
    styleLabel: 'position: absolute; top: -24px; left: 58px',
  },
  {
    id: 'top-right',
    label: 'East',
    styleCheckbox: 'position: absolute; top: 23px; right: -20px',
    styleLabel: 'position: absolute; top: 21px; right: -53px',
  },
  {
    id: 'bottom-left',
    label: 'South',
    styleCheckbox: 'position: absolute; bottom: -20px; left: 35px',
    styleLabel: 'position: absolute; bottom: -23px; left: 56px',
  },
  {
    id: 'bottom-right',
    label: 'West',
    styleCheckbox: 'position: absolute; bottom: 23px; left: -20px',
    styleLabel: 'position: absolute; bottom: 20px; left: -60px',
  },
]

export const csvColumnsNames = ['NorthBarcode', 'SouthBarcode', 'EastBarcode', 'WestBarcode']

export const mapCommandOptions: Options = {
  // probably there is a way to get this data dynamically from the backend
  machine: {
    type: 'str',
    choices: ['echo', 'm1000'],
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
    label: 'For echo mapping only: a yml file with the column headers, otherwise default headers are used',
    required: false,
  },

  experiment_name: {
    type: 'str',
    label: 'Experiment name',
    required: true,
  },
}

export const importCommandOptions: Options = {
  // probably there is a way to get this data dynamically from the backend
  what: {
    type: 'str',
    choices: ['sdf', 'library_plate'],
    label: 'What to import: sdf |  library_plate',
    required: true,
  },
  input_file: {
    type: 'str',
    label: 'The input file',
    required: true,
  },

  library_name: {
    type: 'str',
    label: 'Library name',
    required: false,
  },
  plate_barcode: {
    type: 'str',
    label: 'Plate barcode (for library_plate import only)',
    required: false,
  },
  mapping_file: {
    type: 'str',
    label: 'The mapping file for the sdf columns, otherwise default mapping is used (for sdf import only)',
    required: false,
  },
  // template_name: {
  //   type: 'str',
  //   label: "For template import only: the name of the template, defaults to 'Default'",
  //   required: false,
  // },
}
