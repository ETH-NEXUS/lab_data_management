import { defineStore } from 'pinia'
import { ref } from 'vue'
import type { PaginatedResponse } from '~/types/api'
import type { Plate } from '~/types/lab'
import {
  ADD_CONTROL_LAYOUT_ENDPOINT,
  ADD_CONTROL_LAYOUT_ERROR_MESSAGE,
  CONTROL_PLATES_ENDPOINT,
  CONTROL_PLATES_FETCH_ERROR_MESSAGE,
  CONTROL_PLATES_PARAMS,
  type AddControlLayoutPayload,
} from '~/types/layout'
import { requestApiData, requestApiVoid } from '~/utils/apiRequests'

export const useLayoutStore = defineStore('layoutStore', () => {
  const controlPlates = ref<Plate[]>([])

  const isLoadingControlPlates = ref(false)
  const isAddingControlLayout = ref(false)

  const fetchControlPlates = async (): Promise<Plate[]> => {
    isLoadingControlPlates.value = true
    try {
      const response = await requestApiData<PaginatedResponse<Plate>>(
        CONTROL_PLATES_ENDPOINT,
        {
          method: 'GET',
          params: CONTROL_PLATES_PARAMS,
        },
        CONTROL_PLATES_FETCH_ERROR_MESSAGE,
      )
      controlPlates.value = response.results ?? []
      return controlPlates.value
    } finally {
      isLoadingControlPlates.value = false
    }
  }

  const addControlLayout = async (payload: AddControlLayoutPayload): Promise<void> => {
    isAddingControlLayout.value = true
    try {
      await requestApiVoid(
        ADD_CONTROL_LAYOUT_ENDPOINT,
        {
          method: 'POST',
          body: payload,
        },
        ADD_CONTROL_LAYOUT_ERROR_MESSAGE,
      )
    } finally {
      isAddingControlLayout.value = false
    }
  }

  return {
    controlPlates,
    isLoadingControlPlates,
    isAddingControlLayout,
    fetchControlPlates,
    addControlLayout,
  }
})
