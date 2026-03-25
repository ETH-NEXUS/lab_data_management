import {api} from 'src/boot/axios'
import {defineStore} from 'pinia'
import {harvestProject} from 'src/components/models'

function sortHarvestProjects(harvestProjects: harvestProject[]) {
  const currentYear = new Date().getFullYear().toString()
  harvestProjects.sort((a: harvestProject, b: harvestProject) => {
    const aContainsYear = a.name.includes(currentYear)
    const bContainsYear = b.name.includes(currentYear)
    if (aContainsYear && !bContainsYear) {
      return -1
    } else if (!aContainsYear && bContainsYear) {
      return 1
    } else {
      return a.name.localeCompare(b.name)
    }
  })

  return harvestProjects
}

export const useHarvestStore = defineStore({
  id: 'harvest',
  state: () => ({
    harvestProjects: [] as harvestProject[],
  }),
  actions: {
    async getHarvestProjects() {
      const resp_harvest = await api.get('/api/harvest/harvest_projects/')
      const harvestDataCopy = JSON.parse(JSON.stringify(resp_harvest.data.projects))
      for (let i = 0; i < harvestDataCopy.length; i++) {
        harvestDataCopy[i].value = harvestDataCopy[i].name
        harvestDataCopy[i].label = harvestDataCopy[i].name
      }
      this.harvestProjects = sortHarvestProjects(harvestDataCopy)
    },
    async updateHarvestInfo(projectId: number) {
      const resp = await api.get(`/api/harvest/update_harvest_info/${projectId}/`)
      await this.initialize()
      return resp
    },
    async initialize() {
      if (this.harvestProjects.length === 0) {
        await this.getHarvestProjects()
      }
    },
  },
  persist: {
    enabled: true,
  },
})
