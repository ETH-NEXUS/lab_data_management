import {api} from 'src/boot/axios'
import {FileSystemItem, GeneralFormData} from 'components/models'
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

    async runCommand(formData: GeneralFormData) {
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
    async deleteFile(path: string) {
      await api.post('api/delete_file/', {path: path})
      await this.getDataDirectory()
    },
    async downloadFile(path: string) {
      return await api.post('api/download_file/', {file_path: path})
    },
    async uploadFile(path: string, file: File) {
      const form = new FormData()
      form.append('file', file)
      form.append('directory_path', path)
      const res = await api.post('api/upload_file/', form)
      await this.getDataDirectory()
    },
    getFileContent(path: string) {
      return api.post('api/get_file_content/', {file_path: path})
    },
  },
  persist: {
    enabled: true,
  },
})
