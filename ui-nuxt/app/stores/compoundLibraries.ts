import { defineStore } from 'pinia'
import { ref } from 'vue'
import type { PaginatedResponse } from '~/types/api'
import type { CompoundLibrary, Plate } from '~/types/lab'
import {
  COMPOUND_LIBRARIES_ENDPOINT,
  COMPOUND_LIBRARIES_FETCH_ERROR_MESSAGE,
  COMPOUND_LIBRARY_ADD_PLATE_ENDPOINT,
  COMPOUND_LIBRARY_ADD_PLATE_ERROR_MESSAGE,
  type AddCompoundLibraryPlatePayload,
  type AddCompoundLibraryPlateResult,
} from '~/types/compoundLibraries'
import { requestApiData } from '~/utils/apiRequests'

export const useCompoundLibraryStore = defineStore('compoundLibraryStore', () => {
  const libraries = ref<CompoundLibrary[]>([])
  const isLoadingLibraries = ref(false)
  const isAddingLibraryPlate = ref(false)

  /**
   * Loads compound libraries from backend and stores them in local state.
   *
   * Returned data examples:
   * - `[{ id: 1, name: 'Library A', file_name: 'a.csv', plates: [] }]`
   * - `[]`
   */
  const fetchLibraries = async (): Promise<CompoundLibrary[]> => {
    isLoadingLibraries.value = true
    try {
      const response = await requestApiData<PaginatedResponse<CompoundLibrary>>(
        COMPOUND_LIBRARIES_ENDPOINT,
        {
          method: 'GET',
        },
        COMPOUND_LIBRARIES_FETCH_ERROR_MESSAGE,
      )

      const loadedLibraries: CompoundLibrary[] = []
      const apiLibraries = response.results ?? []

      for (const library of apiLibraries) {
        loadedLibraries.push(library)
      }

      libraries.value = loadedLibraries
      return libraries.value
    } finally {
      isLoadingLibraries.value = false
    }
  }

  /**
   * Creates one new plate and attaches it to one compound library.
   *
   * Accepted input examples:
   * - `libraryId = 5, barcode = 'LIB-001'`
   * - `libraryId = 8, barcode = 'QC-17'`
   */
  const addPlateToLibrary = async (libraryId: number, barcode: string): Promise<AddCompoundLibraryPlateResult> => {
    isAddingLibraryPlate.value = true
    try {
      const payload: AddCompoundLibraryPlatePayload = {
        barcode,
        library: libraryId,
      }

      const createdPlate = await requestApiData<Plate>(
        COMPOUND_LIBRARY_ADD_PLATE_ENDPOINT,
        {
          method: 'POST',
          body: payload,
        },
        COMPOUND_LIBRARY_ADD_PLATE_ERROR_MESSAGE,
      )

      let updatedLibrary: CompoundLibrary | null = null

      for (const library of libraries.value) {
        if (library.id !== libraryId) {
          continue
        }

        const nextPlates: Plate[] = []
        const existingPlates = library.plates ?? []

        for (const plate of existingPlates) {
          nextPlates.push(plate)
        }
        nextPlates.push(createdPlate)

        library.plates = nextPlates
        updatedLibrary = library
        break
      }

      if (!updatedLibrary) {
        await fetchLibraries()

        for (const library of libraries.value) {
          if (library.id === libraryId) {
            updatedLibrary = library
            break
          }
        }
      }

      if (!updatedLibrary) {
        throw new Error(COMPOUND_LIBRARY_ADD_PLATE_ERROR_MESSAGE)
      }

      return {
        updatedLibrary,
        createdPlateId: createdPlate.id,
      }
    } finally {
      isAddingLibraryPlate.value = false
    }
  }

  return {
    libraries,
    isLoadingLibraries,
    isAddingLibraryPlate,
    fetchLibraries,
    addPlateToLibrary,
  }
})
