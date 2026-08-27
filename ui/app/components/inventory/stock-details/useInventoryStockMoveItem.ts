import { computed, ref, watch } from 'vue'
import { useQueryClient } from '@tanstack/vue-query'
import { useInventoryLookupsQuery } from '~/composables/inventory/useInventoryLookupQuery'
import { parsePositiveIntegerList } from '~/components/inventory/add-item/inventorySectorSelection.utils'
import { useInventoryStockStore } from '~/stores/inventory/InventoryStock'
import {
  INVENTORY_STOCK_PAGES_QUERY_KEY,
  INVENTORY_STOCKS_QUERY_KEY,
  type InventoryStockListItem,
} from '~/types/inventory'
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
  const selectedSectorIds = ref<string[]>([])

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
    const nextSectorIds = parsePositiveIntegerList(selectedSectorIds.value).sort((left, right) => left - right)
    const currentSectorIds = (stock?.sectors ?? []).map((sector) => sector.id).sort((left, right) => left - right)

    if (!stock || isSavingMove.value || nextSectorIds.length === 0) {
      return true
    }

    return JSON.stringify(nextSectorIds) === JSON.stringify(currentSectorIds)
  })

  const resetDraft = (): void => {
    const stock = props.stock
    const roomId = stock?.room?.id ?? stock?.sector?.room?.id ?? null
    selectedRoomId.value = roomId ? String(roomId) : ''
    selectedSectorIds.value = (stock?.sectors ?? []).map((sector) => String(sector.id))
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
    const sectorIds = parsePositiveIntegerList(selectedSectorIds.value)

    if (!stock || sectorIds.length === 0) {
      toast.add({ title: 'At least one sector is required', color: 'warning', duration: 2500 })
      return
    }

    isSavingMove.value = true
    try {
      await inventoryStockStore.updateStock(stock.id, { sector_ids: sectorIds })
      await queryClient.invalidateQueries({ queryKey: INVENTORY_STOCKS_QUERY_KEY })
      await queryClient.invalidateQueries({ queryKey: INVENTORY_STOCK_PAGES_QUERY_KEY })
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
    const allowedSectorIds = new Set(filteredSectors.value.map((sector) => String(sector.id)))
    selectedSectorIds.value = selectedSectorIds.value.filter((sectorId) => allowedSectorIds.has(sectorId))
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
    selectedSectorIds,
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
