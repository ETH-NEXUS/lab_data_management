import { defineStore } from 'pinia'
import { ref } from 'vue'
import type { PaginatedResponse } from '~/types/api'
import {
  EXPERIMENT_ADD_MEASUREMENT_ERROR_MESSAGE,
  EXPERIMENT_BARCODE_CREATE_ERROR_MESSAGE,
  EXPERIMENT_BARCODE_DELETE_ERROR_MESSAGE,
  EXPERIMENT_CREATE_ERROR_MESSAGE,
  EXPERIMENT_PLATES_FETCH_ERROR_MESSAGE,
  EXPERIMENT_PREFILL_PLATE_INFO_ERROR_MESSAGE,
  EXPERIMENT_SAVE_PLATE_INFO_ERROR_MESSAGE,
  EXPERIMENT_UPDATE_ERROR_MESSAGE,
  type CreateBarcodeSpecificationPayload,
  type CreateExperimentPayload,
  type Experiment,
  type ExperimentPayload,
  type SavePlateInfoPayload,
} from '~/types/experiments'
import { requestApiData, requestApiVoid } from '~/utils/apiRequests'
import type { CalculatorPayload, Plate, PlateInfo } from '~/types/lab'

export const useExperimentStore = defineStore('experimentStore', () => {
  const isCreatingBarcodeSpecification = ref(false)
  const isDeletingBarcodeSpecification = ref(false)
  const isCreatingExperiment = ref(false)
  const isUpdatingExperiment = ref(false)
  const isLoadingPrefilledPlateInfo = ref(false)
  const isSavingPlateInfo = ref(false)
  const isAddingMeasurement = ref(false)
  const isLoadingExperimentPlates = ref(false)
  const prefilledPlateInfo = ref<PlateInfo[]>([])

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

  /**
   * Loads prefilled experiment plate metadata rows.
   *
   * Accepted input examples:
   * - `experimentId = 21`
   * - `experimentId = 7`
   *
   * Response body example:
   * - `{ plate_info: [{ measurement_label: 'OD600', measurement_timestamp: ['2026-01-01T00:00:00Z'], plate_barcode: 'A001', lib_plate_barcode: 'LIB-1', replicate: '1', cell_type: 'HeLa', condition: 'Control' }] }`
   */
  const fetchPrefilledPlateInfo = async (experimentId: number): Promise<PlateInfo[]> => {
    isLoadingPrefilledPlateInfo.value = true
    try {
      const response = await requestApiData<{ plate_info: PlateInfo[] }>(
        'prefillPlateInfo/',
        {
          method: 'GET',
          params: {
            experiment_id: String(experimentId),
          },
        },
        EXPERIMENT_PREFILL_PLATE_INFO_ERROR_MESSAGE,
      )

      prefilledPlateInfo.value = response.plate_info ?? []
      return prefilledPlateInfo.value
    } finally {
      isLoadingPrefilledPlateInfo.value = false
    }
  }

  /**
   * Saves edited plate metadata rows and refreshes stored rows from backend.
   *
   * Accepted payload example:
   * - `{ experiment_id: 21, plate_info: [{ measurement_label: 'OD600', measurement_timestamp: ['2026-01-01T00:00:00Z'], plate_barcode: 'A001', lib_plate_barcode: 'LIB-1', replicate: '1', cell_type: 'HeLa', condition: 'Control' }] }`
   *
   * Returned value example:
   * - `'success'`
   */
  const savePlateInfo = async (payload: SavePlateInfoPayload): Promise<'success'> => {
    isSavingPlateInfo.value = true
    try {
      await requestApiData<unknown>(
        'save_plate_info/',
        {
          method: 'POST',
          body: payload,
        },
        EXPERIMENT_SAVE_PLATE_INFO_ERROR_MESSAGE,
      )

      await fetchPrefilledPlateInfo(payload.experiment_id)
      return 'success'
    } finally {
      isSavingPlateInfo.value = false
    }
  }

  /**
   * Adds a new derived measurement for one plate or one experiment.
   *
   * Accepted input examples:
   * - `plateId = null`, `experimentId = 21`, `expression = 'A+B'`, `newLabel = 'Combined'`, `usedLabels = ['A', 'B']`
   * - `plateId = 11`, `experimentId = null`, `expression = 'log(A)'`, `newLabel = 'LogA'`, `usedLabels = ['A']`
   *
   * Endpoint example:
   * - `POST /api/plates/null/add_new_measurement/` (experiment calculation)
   * - `POST /api/plates/11/add_new_measurement/` (plate calculation)
   */
  const addNewMeasurement = async (
    plateId: number | null,
    expression: string,
    newLabel: string,
    usedLabels: string[],
    experimentId: number | null = null,
  ): Promise<void> => {
    isAddingMeasurement.value = true
    try {
      const payload: CalculatorPayload = {
        expression,
        new_label: newLabel,
        used_labels: usedLabels,
        separate_time_series_points: false,
      }

      if (plateId !== null) {
        payload.plate_id = plateId
      } else if (experimentId !== null) {
        payload.experiment_id = experimentId
      }

      // Keep old behavior: if a label contains "-->", backend calculates each time series point separately.
      if (usedLabels.some((label) => label.includes('-->'))) {
        payload.separate_time_series_points = true
      }

      await requestApiData<unknown>(
        `plates/${plateId}/add_new_measurement/`,
        {
          method: 'POST',
          body: payload,
        },
        EXPERIMENT_ADD_MEASUREMENT_ERROR_MESSAGE,
      )
    } finally {
      isAddingMeasurement.value = false
    }
  }

  /**
   * Loads all plates belonging to one experiment.
   *
   * Accepted input examples:
   * - `experimentId = 21`
   * - `experimentId = 8`
   *
   * Returned data example:
   * - `[{ id: 11, barcode: 'A001', dimension: { id: 1, name: '384', cols: 24, rows: 16 }, details: { ... }, wells: [] }]`
   */
  const fetchExperimentPlates = async (experimentId: number): Promise<Plate[]> => {
    isLoadingExperimentPlates.value = true
    try {
      const response = await requestApiData<PaginatedResponse<Plate>>(
        'plates/',
        {
          method: 'GET',
          params: {
            experiment: String(experimentId),
          },
        },
        EXPERIMENT_PLATES_FETCH_ERROR_MESSAGE,
      )

      return response.results ?? []
    } finally {
      isLoadingExperimentPlates.value = false
    }
  }

  return {
    isCreatingBarcodeSpecification,
    isDeletingBarcodeSpecification,
    isCreatingExperiment,
    isUpdatingExperiment,
    isLoadingPrefilledPlateInfo,
    isSavingPlateInfo,
    isAddingMeasurement,
    isLoadingExperimentPlates,
    prefilledPlateInfo,
    createBarcodeSpecification,
    deleteBarcodeSpecification,
    createExperiment,
    updateExperiment,
    fetchPrefilledPlateInfo,
    savePlateInfo,
    addNewMeasurement,
    fetchExperimentPlates,
  }
})
