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
  plotView: boolean
  squareCompoundType: boolean
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
  showExperimentResults: boolean
  include_depricated_functionality: boolean
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
        plotView: false,
        squareCompoundType: false,
      },
      navigationTree: {
        filter: '',
      },
      showExperimentResults: false,
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
      include_depricated_functionality: false,
    } as Settings),
  persist: {
    enabled: true,
  },
})
