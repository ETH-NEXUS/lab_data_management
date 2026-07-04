import { computed, ref, watch } from 'vue'
import { useQueryClient } from '@tanstack/vue-query'
import { useExperimentsQuery } from '~/composables/useExperimentQuery'
import { useProjectsQuery } from '~/composables/useProjectsQuery'
import { useInventoryUsageStore } from '~/stores/inventory/InventoryUsageStore'
import {
  INVENTORY_STOCK_PAGES_QUERY_KEY,
  INVENTORY_STOCKS_QUERY_KEY,
  INVENTORY_USAGES_QUERY_KEY,
  type InventoryMaterialDetail,
  type InventoryStockListItem,
} from '~/types/inventory'
import { getErrorMessage } from '~/utils/errors'

type RecordUsageProps = {
  open: boolean
  stock: InventoryStockListItem | null
  materialDetail: InventoryMaterialDetail | null
}

type UsageUnitOption = {
  id: number
  label: string
}

export const useInventoryStockRecordUsage = (props: Readonly<RecordUsageProps>) => {
  const toast = useToast()
  const queryClient = useQueryClient()
  const inventoryUsageStore = useInventoryUsageStore()

  const projectsQuery = useProjectsQuery()
  const experimentsQuery = useExperimentsQuery()

  const isRecordingUsage = ref(false)
  const isSavingUsage = ref(false)
  const selectedProjectId = ref('')
  const selectedExperimentId = ref('')
  const quantityUsed = ref('1')
  const selectedUsageUnitId = ref('')
  const usageNotes = ref('')

  const projects = computed(() => projectsQuery.data.value ?? [])
  const experiments = computed(() => experimentsQuery.data.value ?? [])

  const filteredExperiments = computed(() => {
    const projectId = Number.parseInt(selectedProjectId.value, 10)
    if (!Number.isInteger(projectId) || projectId <= 0) {
      return experiments.value
    }
    return experiments.value.filter((experiment) => experiment.project === projectId)
  })

  const usageUnitOptions = computed<UsageUnitOption[]>(() => {
    const materialUnits = props.materialDetail?.units ?? []
    if (materialUnits.length > 0) {
      return materialUnits.map((unit) => ({
        id: unit.id,
        label: unit.display_name || unit.unit?.label || unit.unit?.name || `Unit #${unit.id}`,
      }))
    }

    const stock = props.stock
    if (!stock) return []

    return [
      {
        id: stock.stock_unit.id,
        label:
          stock.stock_unit.display_name ||
          stock.stock_unit.unit?.label ||
          stock.stock_unit.unit?.name ||
          `Unit #${stock.stock_unit.id}`,
      },
    ]
  })

  const stockItemLabel = computed(() => {
    const stock = props.stock
    if (!stock) return '—'
    return stock.material.label || stock.material.product_name || `Stock #${stock.id}`
  })

  const resetDraft = (): void => {
    selectedProjectId.value = ''
    selectedExperimentId.value = ''
    quantityUsed.value = '1'
    usageNotes.value = ''
    const stockUnitId = props.stock?.stock_unit?.id
    const fallbackUnitId = usageUnitOptions.value[0]?.id ?? null
    selectedUsageUnitId.value = stockUnitId ? String(stockUnitId) : fallbackUnitId ? String(fallbackUnitId) : ''
  }

  const openUsageMode = (): void => {
    if (!props.stock) return
    resetDraft()
    isRecordingUsage.value = true
  }

  const cancelUsageMode = (): void => {
    if (isSavingUsage.value) return
    isRecordingUsage.value = false
    resetDraft()
  }

  const saveUsage = async (): Promise<void> => {
    const stock = props.stock
    if (!stock || isSavingUsage.value) {
      return
    }

    const projectId = Number.parseInt(selectedProjectId.value.trim(), 10)
    if (!Number.isInteger(projectId) || projectId <= 0) {
      toast.add({
        title: 'Project is required',
        color: 'warning',
        duration: 2500,
      })
      return
    }

    const usageUnitId = Number.parseInt(selectedUsageUnitId.value.trim(), 10)
    if (!Number.isInteger(usageUnitId) || usageUnitId <= 0) {
      toast.add({
        title: 'Usage unit is required',
        color: 'warning',
        duration: 2500,
      })
      return
    }

    const quantityText = quantityUsed.value.trim()
    const quantityValue = Number.parseFloat(quantityText)
    if (quantityText === '' || !Number.isFinite(quantityValue) || quantityValue <= 0) {
      toast.add({
        title: 'Enter a valid quantity',
        color: 'warning',
        duration: 2500,
      })
      return
    }

    const experimentText = selectedExperimentId.value.trim()
    const experimentId = experimentText === '' ? null : Number.parseInt(experimentText, 10)
    if (experimentId !== null && (!Number.isInteger(experimentId) || experimentId <= 0)) {
      toast.add({
        title: 'Experiment must be valid',
        color: 'warning',
        duration: 2500,
      })
      return
    }

    isSavingUsage.value = true
    try {
      await inventoryUsageStore.createUsage({
        project_id: projectId,
        experiment_id: experimentId,
        inventory_stock_id: stock.id,
        usage_unit_id: usageUnitId,
        quantity_used: quantityText,
        notes: usageNotes.value.trim() === '' ? null : usageNotes.value.trim(),
      })

      await queryClient.invalidateQueries({ queryKey: INVENTORY_USAGES_QUERY_KEY })
      await queryClient.invalidateQueries({ queryKey: INVENTORY_STOCKS_QUERY_KEY })
      await queryClient.invalidateQueries({ queryKey: INVENTORY_STOCK_PAGES_QUERY_KEY })

      isRecordingUsage.value = false
      resetDraft()

      toast.add({
        title: 'Usage recorded',
        color: 'success',
        duration: 2500,
      })
    } catch (err: unknown) {
      toast.add({
        title: 'Failed to record usage',
        description: getErrorMessage(err),
        color: 'error',
        duration: 4000,
      })
    } finally {
      isSavingUsage.value = false
    }
  }

  watch(selectedProjectId, () => {
    const selectedExperiment = Number.parseInt(selectedExperimentId.value, 10)
    if (!Number.isInteger(selectedExperiment) || selectedExperiment <= 0) {
      selectedExperimentId.value = ''
      return
    }
    const hasMatch = filteredExperiments.value.some((experiment) => experiment.id === selectedExperiment)
    if (!hasMatch) {
      selectedExperimentId.value = ''
    }
  })

  watch(
    () => props.stock,
    () => {
      if (!isRecordingUsage.value) resetDraft()
    },
    { immediate: true },
  )

  watch(
    () => props.open,
    (isOpen) => {
      if (!isOpen) return
      isRecordingUsage.value = false
      resetDraft()
    },
  )

  return {
    isRecordingUsage,
    isSavingUsage,
    selectedProjectId,
    selectedExperimentId,
    quantityUsed,
    selectedUsageUnitId,
    usageNotes,
    projects,
    filteredExperiments,
    usageUnitOptions,
    stockItemLabel,
    openUsageMode,
    cancelUsageMode,
    saveUsage,
  }
}
