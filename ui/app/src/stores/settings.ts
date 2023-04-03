import {WellDetails as WellDetailsType, WellInfo} from './../components/models'
import {Palette} from '../components/helpers'
import {defineStore} from 'pinia'

interface PlatePage {
  splitter: number
  selectedWellInfo: WellInfo | undefined
  wellContent: keyof WellDetailsType
  showHeatmap: boolean
  smallerMapView: boolean
  heatmapPalette: Palette
}

interface NavigationTree {
  expandedNodes: Array<string>
  needsUpdate: boolean
}

interface NavigationTreeFilter {
  filter: string
}

interface WellDetails {
  showStructure: boolean
}

interface Settings {
  platePage: PlatePage
  navigationTree: NavigationTreeFilter
  projectNavigationTree: NavigationTree
  libraryNavigationTree: NavigationTree
  templateNavigationTree: NavigationTree
  wellDetails: WellDetails
}

export const useSettingsStore = defineStore('settings', {
  state: () =>
    ({
      platePage: {
        splitter: 50,
        selectedWellInfo: undefined,
        wellContent: 'hr_position',
        showHeatmap: false,
        smallerMapView: false,
        heatmapPalette: {label: 'OrangeRed', value: {from: '#fff7bc', to: '#993404'}},
      },
      navigationTree: {
        filter: '',
      },
      projectNavigationTree: {
        expandedNodes: [],
        needsUpdate: false,
      },
      libraryNavigationTree: {
        expandedNodes: [],
        needsUpdate: false,
      },
      templateNavigationTree: {
        expandedNodes: [],
        needsUpdate: false,
      },
      wellDetails: {
        showStructure: true,
      },
    } as Settings),
  persist: {
    enabled: true,
  },
})
