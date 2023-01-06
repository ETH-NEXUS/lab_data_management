import {defineStore} from 'pinia'
import {CompoundLibrary} from 'src/components/models'
import {api} from 'src/boot/axios'
import {ref} from 'vue'

export const useCompoundLibraryStore = defineStore('compoundLibrary', () => {
  const libraries = ref<Array<CompoundLibrary>>([])

  const initialize = async () => {
    const resp_p = await api.get('/api/compoundlibraries/')
    libraries.value = resp_p.data.results
  }

  return {libraries, initialize}
})
