import { computed, ref, watch } from 'vue'
import { useQueryClient } from '@tanstack/vue-query'
import { useInventoryStockStore } from '~/stores/inventory/InventoryStock'
import { INVENTORY_STOCKS_QUERY_KEY, type InventoryStockListItem } from '~/types/inventory'
import { getErrorMessage } from '~/utils/errors'

type QuickAdjustProps = {
  open: boolean
  stock: InventoryStockListItem | null
}

export const useInventoryStockQuickAdjust = (props: Readonly<QuickAdjustProps>) => {
  const toast = useToast()
  const queryClient = useQueryClient()
  const inventoryStockStore = useInventoryStockStore()

  const isEditingStock = ref(false)
  const isSavingStockAdjustment = ref(false)
  const editQuantity = ref('')
  const editMinimumQuantity = ref('')
  const editNotes = ref('')

  const stockUnitLabel = computed<string>(() => {
    const stock = props.stock
    if (!stock) return '—'
    return stock.stock_unit.display_name || stock.stock_unit.unit?.label || stock.stock_unit.unit?.name || '—'
  })

  const resetDraft = (): void => {
    const stock = props.stock
    editQuantity.value = stock?.quantity ?? ''
    editMinimumQuantity.value = stock?.minimum_quantity ?? ''
    editNotes.value = stock?.notes ?? ''
  }

  const openEditMode = (): void => {
    if (!props.stock) return
    resetDraft()
    isEditingStock.value = true
  }

  const cancelEditMode = (): void => {
    if (isSavingStockAdjustment.value) return
    isEditingStock.value = false
    resetDraft()
  }

  const saveStockAdjustment = async (): Promise<void> => {
    const stock = props.stock
    if (!stock || isSavingStockAdjustment.value) return

    const quantity = editQuantity.value.trim()
    const minimumQuantity = editMinimumQuantity.value.trim()
    if (quantity === '' || minimumQuantity === '') {
      toast.add({ title: 'Quantity and minimum quantity are required', color: 'warning', duration: 2500 })
      return
    }

    isSavingStockAdjustment.value = true
    try {
      await inventoryStockStore.updateStock(stock.id, {
        quantity,
        minimum_quantity: minimumQuantity,
        notes: editNotes.value.trim() === '' ? null : editNotes.value.trim(),
      })
      await queryClient.invalidateQueries({ queryKey: INVENTORY_STOCKS_QUERY_KEY })
      isEditingStock.value = false
      toast.add({ title: 'Stock adjusted', color: 'success', duration: 2500 })
    } catch (err: unknown) {
      toast.add({
        title: 'Failed to adjust stock',
        description: getErrorMessage(err),
        color: 'error',
        duration: 4000,
      })
    } finally {
      isSavingStockAdjustment.value = false
    }
  }

  watch(
    () => props.stock,
    () => {
      if (!isEditingStock.value) resetDraft()
    },
    { immediate: true },
  )

  watch(
    () => props.open,
    (isOpen) => {
      if (!isOpen) return
      isEditingStock.value = false
      resetDraft()
    },
  )

  return {
    isEditingStock,
    isSavingStockAdjustment,
    editQuantity,
    editMinimumQuantity,
    editNotes,
    stockUnitLabel,
    openEditMode,
    cancelEditMode,
    saveStockAdjustment,
  }
}
