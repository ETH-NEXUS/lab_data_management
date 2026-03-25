import { defineStore } from 'pinia'
import { ref } from 'vue'
import type { PaginatedResponse } from '~/types/api'
import type { Plate, PlateLabelValue } from '~/types/lab'
import {
  getMultiplePlatesFoundMessage,
  getPlateNotFoundMessage,
  PLATE_BARCODES_ENDPOINT,
  PLATE_FETCH_ERROR_MESSAGE,
  PLATES_ENDPOINT,
  PLATE_TARGET_BARCODES_FETCH_ERROR_MESSAGE,
  PLATE_TEMPLATE_BARCODES_FETCH_ERROR_MESSAGE,
  type PlateBarcode,
  type PlatePageData,
} from '~/types/plates'
import { requestApiData } from '~/utils/apiRequests'
import { getErrorMessage } from '~/utils/errors'

const normalizePlateBarcode = (barcode: PlateBarcode): string => barcode.trim()

/**
 * Filters out the currently opened plate from a list of selectable plate options.
 *
 * Accepted data example:
 * - options: `[{ label: 'A001', value: 1 }, { label: 'A002', value: 2 }]`
 * - excludedBarcode: `'A001'`
 *
 * Returned data example:
 * - `[{ label: 'A002', value: 2 }]`
 */
const excludeCurrentPlate = (options: PlateLabelValue[], excludedBarcode: string | null): PlateLabelValue[] => {
  if (!excludedBarcode) return [...options]
  return options.filter((option) => option.label !== excludedBarcode)
}

export const usePlateStore = defineStore('plateStore', () => {
  const currentPlateBarcode = ref<string | null>(null)
  const currentPlate = ref<Plate | null>(null)
  const targetPlateBarcodeOptions = ref<PlateLabelValue[]>([])
  const templatePlateBarcodeOptions = ref<PlateLabelValue[]>([])

  const isInitializingPlatePage = ref(false)
  const isLoadingPlate = ref(false)
  const isLoadingTargetPlateBarcodes = ref(false)
  const isLoadingTemplatePlateBarcodes = ref(false)
  const error = ref<string | null>(null)

  /**
   * Clears current plate-page data while keeping store actions available.
   */
  const clearPlatePageState = (): void => {
    currentPlateBarcode.value = null
    currentPlate.value = null
    targetPlateBarcodeOptions.value = []
    templatePlateBarcodeOptions.value = []
    error.value = null
  }

  /**
   * Loads one plate strictly by barcode (legacy behavior).
   *
   * Accepted barcode examples:
   * - `'demo_1'`
   * - `'250401control'`
   *
   * Returned data examples:
   * - `{ id: 42, barcode: 'demo_1', dimension: { id: 1, name: '384', cols: 24, rows: 16 } }`
   * - `null` (not found)
   */
  const fetchPlateByBarcode = async (barcode: PlateBarcode): Promise<Plate | null> => {
    isLoadingPlate.value = true
    const normalizedBarcode = normalizePlateBarcode(barcode)
    currentPlateBarcode.value = normalizedBarcode
    try {
      const data = await requestApiData<PaginatedResponse<Plate>>(
        PLATES_ENDPOINT,
        {
          method: 'GET',
          params: { barcode: normalizedBarcode },
        },
        PLATE_FETCH_ERROR_MESSAGE,
      )

      const results = data.results ?? []

      if (results.length > 1) {
        throw new Error(getMultiplePlatesFoundMessage(normalizedBarcode))
      }

      if (results.length === 1) {
        const firstPlate = results[0]
        if (!firstPlate) {
          currentPlate.value = null
          error.value = getPlateNotFoundMessage(normalizedBarcode)
          return null
        }

        currentPlate.value = firstPlate
        return firstPlate
      }

      currentPlate.value = null
      error.value = getPlateNotFoundMessage(normalizedBarcode)
      return null
    } catch (err: unknown) {
      error.value = getErrorMessage(err)
      throw err
    } finally {
      isLoadingPlate.value = false
    }
  }

  /**
   * Loads selectable target plate barcode options (`experiment=true&library=true`).
   *
   * Accepted argument examples:
   * - `null`
   * - `'demo_1'`
   *
   * Returned data example:
   * - `[{ label: 'A001', value: 11, experiment: { id: 5, name: 'Dose' }, library: null }]`
   */
  const fetchTargetPlateBarcodeOptions = async (currentPlateBarcode: string | null): Promise<PlateLabelValue[]> => {
    isLoadingTargetPlateBarcodes.value = true

    try {
      const data = await requestApiData<PlateLabelValue[]>(
        PLATE_BARCODES_ENDPOINT,
        {
          method: 'GET',
          params: {
            experiment: 'true',
            library: 'true',
          },
        },
        PLATE_TARGET_BARCODES_FETCH_ERROR_MESSAGE,
      )

      targetPlateBarcodeOptions.value = excludeCurrentPlate(data, currentPlateBarcode)
      return targetPlateBarcodeOptions.value
    } catch (err: unknown) {
      error.value = getErrorMessage(err)
      throw err
    } finally {
      isLoadingTargetPlateBarcodes.value = false
    }
  }

  /**
   * Loads selectable template barcode options for the current plate barcode.
   *
   * Accepted argument examples:
   * - `'demo_1'`
   * - `null`
   *
   * Returned data example:
   * - `[{ label: '__TEMPL__Project_Experiment', value: 77, experiment: null, library: null }]`
   */
  const fetchTemplatePlateBarcodeOptions = async (currentPlateBarcode: string | null): Promise<PlateLabelValue[]> => {
    isLoadingTemplatePlateBarcodes.value = true

    try {
      if (!currentPlateBarcode) {
        templatePlateBarcodeOptions.value = []
        return templatePlateBarcodeOptions.value
      }

      const data = await requestApiData<PlateLabelValue[]>(
        PLATE_BARCODES_ENDPOINT,
        {
          method: 'GET',
          params: { barcode: currentPlateBarcode },
        },
        PLATE_TEMPLATE_BARCODES_FETCH_ERROR_MESSAGE,
      )

      templatePlateBarcodeOptions.value = excludeCurrentPlate(data, currentPlateBarcode)
      return templatePlateBarcodeOptions.value
    } catch (err: unknown) {
      error.value = getErrorMessage(err)
      throw err
    } finally {
      isLoadingTemplatePlateBarcodes.value = false
    }
  }

  /**
   * Initializes all plate-page fetch state for a single opened route.
   *
   * Behavior:
   * - loads current plate,
   * - loads both barcode options lists.
   *
   * Accepted barcode examples:
   * - `'demo_1'`
   * - `'250401control'`
   *
   * Returned data example:
   * - `{ plate: { id: 42, barcode: 'demo_1' }, targetPlateBarcodeOptions: [], templatePlateBarcodeOptions: [] }`
   */
  const initializePlatePage = async (barcode: PlateBarcode): Promise<PlatePageData> => {
    isInitializingPlatePage.value = true
    error.value = null

    try {
      const plate = await fetchPlateByBarcode(barcode)
      const currentPlateBarcode = plate?.barcode ?? null

      await Promise.all([
        fetchTargetPlateBarcodeOptions(currentPlateBarcode),
        fetchTemplatePlateBarcodeOptions(currentPlateBarcode),
      ])

      return {
        plate: currentPlate.value,
        targetPlateBarcodeOptions: targetPlateBarcodeOptions.value,
        templatePlateBarcodeOptions: templatePlateBarcodeOptions.value,
      }
    } catch (err: unknown) {
      // Keep thrown behavior so caller can decide UX flow.
      error.value = getErrorMessage(err)
      throw err
    } finally {
      isInitializingPlatePage.value = false
    }
  }

  return {
    currentPlateBarcode,
    currentPlate,
    targetPlateBarcodeOptions,
    templatePlateBarcodeOptions,
    isInitializingPlatePage,
    isLoadingPlate,
    isLoadingTargetPlateBarcodes,
    isLoadingTemplatePlateBarcodes,
    error,
    clearPlatePageState,
    fetchPlateByBarcode,
    fetchTargetPlateBarcodeOptions,
    fetchTemplatePlateBarcodeOptions,
    initializePlatePage,
  }
})
