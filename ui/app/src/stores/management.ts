import {api} from 'src/boot/axios'
import {FileSystemItem, FormData} from 'components/models'
import {defineStore} from 'pinia'

export const useManagementStore = defineStore({
  id: 'management',
  state: () => ({
    dataDirectory: {type: 'directory', name: '', children: [], path: ''} as FileSystemItem,
    selectedPath: '',
    selectedPaths: [] as string[],
    commandOutput: '' as string,
  }),
  actions: {
    async getDataDirectory() {
      const res = await api.get('/api/directory_content/')
      this.dataDirectory = res.data.directory_content
    },

    async initialize() {
      await this.getDataDirectory()
    },

    async runCommand(formData: FormData) {
      const res = await api.post('/api/run_command/', {form_data: formData})
      this.commandOutput = res.data.command_output
    },
  },
  persist: {
    enabled: true,
  },
})
