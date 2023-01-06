import {WellInfo} from './../components/models'
import {defineStore} from 'pinia'

interface PlatePage {
  splitter: number
  selectedWellInfo: WellInfo | undefined
}

interface NavigationTree {
  expandedNodes: Array<string>
  needsUpdate: boolean
}

interface WellDetails {
  showStructure: boolean
}

interface Settings {
  platePage: PlatePage
  navigationTree: NavigationTree
  wellDetails: WellDetails
}

export const useSettingsStore = defineStore('settings', {
  state: () =>
    ({
      platePage: {
        splitter: 50,
        selectedWellInfo: undefined,
      },
      navigationTree: {
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
