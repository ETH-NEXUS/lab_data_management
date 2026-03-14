import { defineStore } from 'pinia'
import { ref } from 'vue'
import { useAPI } from '~/composables/useAPI'
import {
  EXPERIMENT_BARCODE_CREATE_ERROR_MESSAGE,
  EXPERIMENT_CREATE_ERROR_MESSAGE,
  EXPERIMENT_UPDATE_ERROR_MESSAGE,
  type CreateBarcodeSpecificationPayload,
  type CreateExperimentPayload,
  type Experiment,
  type ExperimentPayload,
} from '~/types/experiments'

export const useExperimentStore = defineStore('experimentStore', () => {
  const isCreatingBarcodeSpecification = ref(false)
  const isCreatingExperiment = ref(false)
  const isUpdatingExperiment = ref(false)

  /**
   * Creates a new experiment in the backend.
   *
   * Accepted payload examples:
   * - `{ name: 'Dose response', project: 5 }`
   * - `{ name: 'QC run', project: 8 }`
   *
   * Returned data example:
   * - `{ id: 21, name: 'Dose response', project: 5 }`
   */
  const createExperiment = async (payload: CreateExperimentPayload): Promise<Experiment> => {
    isCreatingExperiment.value = true
    try {
      const { data, error } = await useAPI<Experiment>('experiments/', {
        method: 'POST',
        body: payload,
      })

      if (error.value || !data.value) {
        throw (error.value ?? new Error(EXPERIMENT_CREATE_ERROR_MESSAGE)) as Error
      }

      return data.value
    } finally {
      isCreatingExperiment.value = false
    }
  }

  /**
   * Creates a barcode specification for one experiment.
   *
   * Accepted payload examples:
   * - `{ experiment_id: 21, prefix: 'A001', number_of_plates: 12, sides: ['North', 'South'] }`
   * - `{ experiment_id: 22, prefix: 'B', number_of_plates: 4, sides: ['East'] }`
   *
   * Backend response example:
   * - HTTP 200 with an empty body (no JSON payload)
   */
  const createBarcodeSpecification = async (payload: CreateBarcodeSpecificationPayload): Promise<void> => {
    isCreatingBarcodeSpecification.value = true
    try {
      const { error } = await useAPI<unknown>('experiments/barcodes/', {
        method: 'POST',
        body: payload,
      })

      // This endpoint is successful even when response body is empty.
      if (error.value) {
        throw (error.value ?? new Error(EXPERIMENT_BARCODE_CREATE_ERROR_MESSAGE)) as Error
      }
    } finally {
      isCreatingBarcodeSpecification.value = false
    }
  }

  /**
   * Updates an experiment on the backend.
   *
   * Accepted payload example:
   * - `{ name: 'Updated experiment name' }`
   *
   * Returned data example:
   * - `{ id: 21, name: 'Updated experiment name', project: 5 }`
   */
  const updateExperiment = async (experimentId: number, payload: ExperimentPayload): Promise<Experiment> => {
    isUpdatingExperiment.value = true
    try {
      const { data, error } = await useAPI<Experiment>(`experiments/${experimentId}/`, {
        method: 'PATCH',
        body: payload,
      })

      if (error.value || !data.value) {
        throw (error.value ?? new Error(EXPERIMENT_UPDATE_ERROR_MESSAGE)) as Error
      }

      return data.value
    } finally {
      isUpdatingExperiment.value = false
    }
  }

  return {
    isCreatingBarcodeSpecification,
    isCreatingExperiment,
    isUpdatingExperiment,
    createBarcodeSpecification,
    createExperiment,
    updateExperiment,
  }
})
