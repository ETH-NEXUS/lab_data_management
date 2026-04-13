import { computed, ref, watch } from 'vue'
import { useQueryClient } from '@tanstack/vue-query'
import { useInventoryLookupsQuery } from '~/composables/inventory/useInventoryLookupQuery'
import { useInventoryStockStore } from '~/stores/inventory/InventoryStock'
import { INVENTORY_STOCKS_QUERY_KEY, type InventoryStockListItem } from '~/types/inventory'
import { getErrorMessage } from '~/utils/errors'

type MoveItemProps = {
  open: boolean
  stock: InventoryStockListItem | null
}

export const useInventoryStockMoveItem = (props: Readonly<MoveItemProps>) => {
  const toast = useToast()
  const queryClient = useQueryClient()
  const inventoryStockStore = useInventoryStockStore()
  const lookupsQuery = useInventoryLookupsQuery()

  const isMovingStock = ref(false)
  const isSavingMove = ref(false)
  const selectedRoomId = ref('')
  const selectedSectorId = ref('')

  const rooms = computed(() => lookupsQuery.data.value?.rooms ?? [])
  const filteredSectors = computed(() => {
    const roomId = Number.parseInt(selectedRoomId.value, 10)
    if (!Number.isInteger(roomId) || roomId <= 0) return []
    return (lookupsQuery.data.value?.sectors ?? []).filter((sector) => sector.room?.id === roomId)
  })
  const isLookupsLoading = computed<boolean>(() => lookupsQuery.isPending.value)
  const lookupsErrorMessage = computed<string | null>(() => {
    const err = lookupsQuery.error.value
    if (!err) return null
    return getErrorMessage(err)
  })
  const isMoveSaveDisabled = computed<boolean>(() => {
    const stock = props.stock
    if (!stock || isSavingMove.value || selectedSectorId.value.trim() === '') return true
    return String(stock.sector.id) === selectedSectorId.value
  })

  const resetDraft = (): void => {
    const stock = props.stock
    const roomId = stock?.room?.id ?? stock?.sector?.room?.id ?? null
    selectedRoomId.value = roomId ? String(roomId) : ''
    selectedSectorId.value = stock?.sector?.id ? String(stock.sector.id) : ''
  }

  const openMoveMode = (): void => {
    if (!props.stock) return
    resetDraft()
    isMovingStock.value = true
  }

  const cancelMoveMode = (): void => {
    if (isSavingMove.value) return
    isMovingStock.value = false
    resetDraft()
  }

  const saveMove = async (): Promise<void> => {
    const stock = props.stock
    const sectorId = Number.parseInt(selectedSectorId.value, 10)
    if (!stock || !Number.isInteger(sectorId) || sectorId <= 0) {
      toast.add({ title: 'Sector is required', color: 'warning', duration: 2500 })
      return
    }

    isSavingMove.value = true
    try {
      await inventoryStockStore.updateStock(stock.id, { sector_id: sectorId })
      await queryClient.invalidateQueries({ queryKey: INVENTORY_STOCKS_QUERY_KEY })
      isMovingStock.value = false
      toast.add({ title: 'Location updated', color: 'success', duration: 2500 })
    } catch (err: unknown) {
      toast.add({
        title: 'Failed to move item',
        description: getErrorMessage(err),
        color: 'error',
        duration: 4000,
      })
    } finally {
      isSavingMove.value = false
    }
  }

  watch(selectedRoomId, () => {
    const currentSectorId = Number.parseInt(selectedSectorId.value, 10)
    if (!Number.isInteger(currentSectorId) || currentSectorId <= 0) {
      selectedSectorId.value = ''
      return
    }
    const hasCurrentSector = filteredSectors.value.some((sector) => sector.id === currentSectorId)
    if (!hasCurrentSector) selectedSectorId.value = ''
  })

  watch(
    () => props.stock,
    () => !isMovingStock.value && resetDraft(),
    { immediate: true },
  )
  watch(
    () => props.open,
    (isOpen) => isOpen && cancelMoveMode(),
  )

  return {
    isMovingStock,
    isSavingMove,
    selectedRoomId,
    selectedSectorId,
    rooms,
    filteredSectors,
    isLookupsLoading,
    lookupsErrorMessage,
    isMoveSaveDisabled,
    openMoveMode,
    cancelMoveMode,
    saveMove,
  }
}
