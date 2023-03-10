import {WellInfo} from './../components/models'
import {defineStore} from 'pinia'
import {Palette} from 'components/data'
import {ref} from 'vue'

interface PlatePage {
  splitter: number
  selectedWellInfo: WellInfo | undefined
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
  palette: Palette
  showHeatmap: boolean
}

export const useSettingsStore = defineStore('settings', {
  state: () =>
    ({
      platePage: {
        splitter: 50,
        selectedWellInfo: undefined,
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
      palette: {value: 'orange', label: 'Palette 1', from: '#fff7bc', to: '#993404'},

      showHeatmap: false,
    } as Settings),
  persist: {
    enabled: true,
  },
})
