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
      this.startLongPolling(formData['room_name'].toString()) // no await here
      await api.post('/api/run_command/', {form_data: formData})
    },

    async startLongPolling(roomName: string) {
      try {
        const res = await api.get(`api/long_polling/${roomName}/`)
        const message = res.data.message
        const status = res.data.status

        if (message) {
          this.commandOutput += '\n' + message
        }

        if (status !== 'completed') {
          setTimeout(() => this.startLongPolling(roomName), 500)
        }
      } catch (error) {
        console.error('Long polling error:', error)
      }
    },
  },
  persist: {
    enabled: true,
  },
})
