import { defineStore } from 'pinia'
import { ref } from 'vue'
import {
  EXPERIMENT_BARCODE_CREATE_ERROR_MESSAGE,
  EXPERIMENT_BARCODE_DELETE_ERROR_MESSAGE,
  EXPERIMENT_CREATE_ERROR_MESSAGE,
  EXPERIMENT_UPDATE_ERROR_MESSAGE,
  type CreateBarcodeSpecificationPayload,
  type CreateExperimentPayload,
  type Experiment,
  type ExperimentPayload,
} from '~/types/experiments'
import { requestApiData, requestApiVoid } from '~/utils/apiRequests'

export const useExperimentStore = defineStore('experimentStore', () => {
  const isCreatingBarcodeSpecification = ref(false)
  const isDeletingBarcodeSpecification = ref(false)
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
      const experiment = await requestApiData<Experiment>(
        'experiments/',
        {
          method: 'POST',
          body: payload,
        },
        EXPERIMENT_CREATE_ERROR_MESSAGE,
      )

      return experiment
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
      await requestApiVoid(
        'experiments/barcodes/',
        {
          method: 'POST',
          body: payload,
        },
        EXPERIMENT_BARCODE_CREATE_ERROR_MESSAGE,
      )
    } finally {
      isCreatingBarcodeSpecification.value = false
    }
  }

  /**
   * Deletes one barcode specification by id.
   *
   * Accepted id examples:
   * - `7`
   * - `42`
   *
   * Endpoint example:
   * - `DELETE /api/barcodespecifications/7/`
   */
  const deleteBarcodeSpecification = async (barcodeSpecificationId: number): Promise<void> => {
    isDeletingBarcodeSpecification.value = true
    try {
      await requestApiVoid(
        `barcodespecifications/${barcodeSpecificationId}/`,
        { method: 'DELETE' },
        EXPERIMENT_BARCODE_DELETE_ERROR_MESSAGE,
      )
    } finally {
      isDeletingBarcodeSpecification.value = false
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
      const experiment = await requestApiData<Experiment>(
        `experiments/${experimentId}/`,
        {
          method: 'PATCH',
          body: payload,
        },
        EXPERIMENT_UPDATE_ERROR_MESSAGE,
      )

      return experiment
    } finally {
      isUpdatingExperiment.value = false
    }
  }

  return {
    isCreatingBarcodeSpecification,
    isDeletingBarcodeSpecification,
    isCreatingExperiment,
    isUpdatingExperiment,
    createBarcodeSpecification,
    deleteBarcodeSpecification,
    createExperiment,
    updateExperiment,
  }
})
