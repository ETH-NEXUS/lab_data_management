import { computed, ref, watch } from 'vue'
import { useExperimentsQuery } from '~/composables/useExperimentQuery'
import { useProjectsQuery } from '~/composables/useProjectsQuery'
import type { InventoryMaterialDetail, InventoryStockListItem } from '~/types/inventory'

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

  const saveUsage = async (): Promise<void> => {}

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
