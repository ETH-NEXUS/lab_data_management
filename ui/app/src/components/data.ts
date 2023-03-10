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

type Palette = {
  val: string
  label: string
  from: string
  to: string
}

export const palettes: Record<string, Palette> = {
  green_red: {val: 'green_red', label: 'Palette 1', from: '#00FF00', to: '#FF0000'},
  green: {val: 'green', label: 'Palette 2', from: '#c7e9c0', to: '#006d2c'},
  blue: {val: 'blue', label: 'Palette 3', from: '#92c5de', to: '#0b2746'},
  orange: {val: 'orange', label: 'Palette 4', from: '#fff7bc', to: '#993404'},
  green_brown: {val: 'green_brown', label: 'Palette 5', from: '#b8e186', to: '#662506'},
}
