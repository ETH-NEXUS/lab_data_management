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

  const addPlate = async (library: CompoundLibrary, barcode: string) => {
    const resp = await api.post('/api/plates/', {
      barcode: barcode,
      library: library.id,
    })
    const plate = resp.data
    library.plates.push(plate)
    return plate
  }

  return {libraries, initialize, addPlate}
})
