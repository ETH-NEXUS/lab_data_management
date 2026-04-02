<script setup lang="ts">
import { computed, onBeforeUnmount, onMounted, ref, watch } from 'vue'
import type {
  InventoryMaterialDetail,
  InventoryOrderListItem,
  InventoryStockListItem,
  InventoryUsageListItem,
} from '~/types/inventory'
import { formatDateTime } from '~/utils/dateTime'

type DrawerField = {
  label: string
  value: string
  wide?: boolean
}

type Props = {
  open: boolean
  stock: InventoryStockListItem | null
  materialDetail: InventoryMaterialDetail | null
  usageEntries: InventoryUsageListItem[]
  orderEntries: InventoryOrderListItem[]
  isMaterialLoading?: boolean
  isUsagesLoading?: boolean
  isOrdersLoading?: boolean
}

const props = withDefaults(defineProps<Props>(), {
  isMaterialLoading: false,
  isUsagesLoading: false,
  isOrdersLoading: false,
})

const emit = defineEmits<{
  (e: 'close'): void
}>()

const { t } = useI18n()

const isMetadataExpanded = ref(false)

watch(
  () => props.open,
  (isOpen) => {
    if (isOpen) {
      isMetadataExpanded.value = false
    }
  },
)

const onEscape = (event: KeyboardEvent): void => {
  if (event.key === 'Escape' && props.open) {
    emit('close')
  }
}

onMounted(() => {
  window.addEventListener('keydown', onEscape)
})

onBeforeUnmount(() => {
  window.removeEventListener('keydown', onEscape)
})

/**
 * Converts optional values into a readable fallback-safe string.
 *
 * Input examples:
 * - `'Vendor A'`
 * - `null`
 * - `''`
 *
 * Returned examples:
 * - `'Vendor A'`
 * - `'—'`
 */
const toDisplayValue = (value: unknown): string => {
  if (value == null) {
    return t('inventory.stock_drawer.values.none')
  }

  const textValue = String(value).trim()
  if (textValue === '') {
    return t('inventory.stock_drawer.values.none')
  }

  return textValue
}

/**
 * Resolves one room/sector location label from stock values.
 *
 * Accepted stock example:
 * - `{ room: { label: 'Room A' }, sector: { label: 'Shelf 3' } }`
 *
 * Returned examples:
 * - `'Room A / Shelf 3'`
 * - `'—'`
 */
const locationLabel = computed<string>(() => {
  const stock = props.stock
  if (!stock) {
    return t('inventory.stock_drawer.values.none')
  }

  const roomLabel = stock.room?.label ?? stock.room?.name ?? stock.sector?.room?.label ?? stock.sector?.room?.name ?? ''
  const sectorLabel = stock.sector?.label ?? stock.sector?.name ?? ''

  if (roomLabel !== '' && sectorLabel !== '') {
    return `${roomLabel} / ${sectorLabel}`
  }

  if (roomLabel !== '') {
    return roomLabel
  }

  if (sectorLabel !== '') {
    return sectorLabel
  }

  return t('inventory.stock_drawer.values.none')
})

const inventoryStatusLabel = computed<string>(() => {
  if (!props.stock) {
    return t('inventory.stock_drawer.values.none')
  }
  return props.stock.inventory_status === 'low'
    ? t('inventory.stock_table.status_labels.low')
    : t('inventory.stock_table.status_labels.in_stock')
})

const inventoryStatusColor = computed<'warning' | 'success'>(() => {
  return props.stock?.inventory_status === 'low' ? 'warning' : 'success'
})

const stockUnitLabel = computed<string>(() => {
  if (!props.stock) {
    return t('inventory.stock_drawer.values.none')
  }
  return toDisplayValue(
    props.stock.stock_unit.display_name || props.stock.stock_unit.unit?.label || props.stock.stock_unit.unit?.name,
  )
})

const quantityLabel = computed<string>(() => {
  if (!props.stock) {
    return t('inventory.stock_drawer.values.none')
  }
  return `${toDisplayValue(props.stock.quantity)} ${stockUnitLabel.value}`.trim()
})

const minimumQuantityLabel = computed<string>(() => {
  if (!props.stock) {
    return t('inventory.stock_drawer.values.none')
  }
  return `${toDisplayValue(props.stock.minimum_quantity)} ${stockUnitLabel.value}`.trim()
})

const baseUnitLabel = computed<string>(() => {
  const stock = props.stock
  if (!stock) {
    return t('inventory.stock_drawer.values.none')
  }

  // Use material detail units when available.
  const materialUnits = props.materialDetail?.units ?? []
  const baseUnit = materialUnits.find((unit) => unit.is_base_unit)
  if (baseUnit) {
    return toDisplayValue(baseUnit.display_name || baseUnit.unit?.label || baseUnit.unit?.name)
  }

  // Fallback: stock unit can already be base.
  if (stock.stock_unit.is_base_unit) {
    return toDisplayValue(stock.stock_unit.display_name || stock.stock_unit.unit?.label || stock.stock_unit.unit?.name)
  }

  // Fallback to stock unit label.
  return toDisplayValue(stock.stock_unit.unit?.label || stock.stock_unit.unit?.name)
})

const attributeLabel = computed<string>(() => {
  const stock = props.stock
  if (!stock) {
    return t('inventory.stock_drawer.values.none')
  }

  const labels: string[] = []
  for (const attribute of stock.material.attributes ?? []) {
    const label = attribute.label || attribute.name
    if (label && label.trim() !== '') {
      labels.push(label)
    }
  }

  if (labels.length === 0) {
    return t('inventory.stock_drawer.values.no_attributes')
  }

  return labels.join(', ')
})

const operationalFields = computed<DrawerField[]>(() => {
  const stock = props.stock
  if (!stock) return []

  return [
    {
      label: t('inventory.stock_drawer.fields.product_name'),
      value: toDisplayValue(stock.material.product_name),
      wide: true,
    },
    {
      label: t('inventory.stock_drawer.fields.quantity'),
      value: quantityLabel.value,
    },
    {
      label: t('inventory.stock_drawer.fields.minimum_quantity'),
      value: minimumQuantityLabel.value,
    },
    {
      label: t('inventory.stock_drawer.fields.location'),
      value: locationLabel.value,
      wide: true,
    },
    {
      label: t('inventory.stock_drawer.fields.lot_number'),
      value: toDisplayValue(stock.lot_number),
    },
    {
      label: t('inventory.stock_drawer.fields.expiry_date'),
      value: stock.expiry_date
        ? formatDateTime(stock.expiry_date, { dateStyle: 'medium' }, t('inventory.stock_drawer.values.none'))
        : t('inventory.stock_drawer.values.none'),
    },
    {
      label: t('inventory.stock_drawer.fields.notes'),
      value: stock.notes || t('inventory.stock_drawer.values.no_notes'),
      wide: true,
    },
  ]
})

const unitConversionFields = computed<DrawerField[]>(() => {
  const stock = props.stock
  if (!stock) return []

  return [
    {
      label: t('inventory.stock_drawer.fields.stock_unit'),
      value: stockUnitLabel.value,
    },
    {
      label: t('inventory.stock_drawer.fields.base_unit'),
      value: baseUnitLabel.value,
    },
    {
      label: t('inventory.stock_drawer.fields.base_units_per_unit'),
      value: toDisplayValue(stock.stock_unit.base_units_per_unit),
    },
  ]
})

const materialIdentityFields = computed<DrawerField[]>(() => {
  const stock = props.stock
  if (!stock) return []

  return [
    {
      label: t('inventory.stock_drawer.fields.device_type'),
      value: toDisplayValue(stock.material.device_type?.label || stock.material.device_type?.name),
    },
    {
      label: t('inventory.stock_drawer.fields.item_type'),
      value: toDisplayValue(stock.material.item_type?.label || stock.material.item_type?.name),
    },
    {
      label: t('inventory.stock_drawer.fields.attributes'),
      value: attributeLabel.value,
      wide: true,
    },
    {
      label: t('inventory.stock_drawer.fields.capacity_value'),
      value: toDisplayValue(stock.material.capacity_value),
    },
    {
      label: t('inventory.stock_drawer.fields.capacity_unit'),
      value: toDisplayValue(stock.material.capacity_unit),
    },
    {
      label: t('inventory.stock_drawer.fields.description'),
      value: toDisplayValue(props.materialDetail?.description),
      wide: true,
    },
  ]
})

const supplierFields = computed<DrawerField[]>(() => {
  const stock = props.stock
  if (!stock) return []

  return [
    {
      label: t('inventory.stock_drawer.fields.manufacturer'),
      value: toDisplayValue(stock.material.manufacturer?.label || stock.material.manufacturer?.name),
    },
    {
      label: t('inventory.stock_drawer.fields.vendor'),
      value: toDisplayValue(stock.material.vendor?.label || stock.material.vendor?.name),
    },
    {
      label: t('inventory.stock_drawer.fields.manufacturer_catalog_number'),
      value: toDisplayValue(stock.material.manufacturer_catalog_number),
      wide: true,
    },
    {
      label: t('inventory.stock_drawer.fields.vendor_catalog_number'),
      value: toDisplayValue(stock.material.vendor_catalog_number),
      wide: true,
    },
    {
      label: t('inventory.stock_drawer.fields.brand'),
      value: toDisplayValue(stock.material.brand?.label || stock.material.brand?.name),
    },
    {
      label: t('inventory.stock_drawer.fields.default_cost'),
      value: toDisplayValue(stock.material.default_cost),
    },
  ]
})

const sortedUsageEntries = computed<InventoryUsageListItem[]>(() => {
  const usages = [...props.usageEntries]
  usages.sort((leftEntry, rightEntry) => {
    const leftDate = new Date(leftEntry.used_at).getTime()
    const rightDate = new Date(rightEntry.used_at).getTime()
    return rightDate - leftDate
  })
  return usages
})

const linkedExperimentProjectLabels = computed<string[]>(() => {
  const labels = new Set<string>()

  for (const usage of props.usageEntries) {
    const projectLabel = usage.project?.label || usage.project?.name || ''
    const experimentLabel = usage.experiment?.label || usage.experiment?.name || ''

    if (projectLabel !== '' && experimentLabel !== '') {
      labels.add(`${projectLabel} / ${experimentLabel}`)
      continue
    }

    if (projectLabel !== '') {
      labels.add(projectLabel)
      continue
    }

    if (experimentLabel !== '') {
      labels.add(experimentLabel)
    }
  }

  return Array.from(labels)
})

const sortedOrderEntries = computed<InventoryOrderListItem[]>(() => {
  const orders = [...props.orderEntries]
  orders.sort((leftOrder, rightOrder) => {
    const leftDate = new Date(leftOrder.order_date).getTime()
    const rightDate = new Date(rightOrder.order_date).getTime()
    return rightDate - leftDate
  })
  return orders
})

const responsibleUsers = computed<string[]>(() => {
  const labels = new Set<string>()

  for (const order of props.orderEntries) {
    const userLabel = order.who_ordered?.label || order.who_ordered?.full_name || order.who_ordered?.username || ''
    if (userLabel !== '') {
      labels.add(userLabel)
    }
  }

  return Array.from(labels)
})

const metadataFields = computed<DrawerField[]>(() => {
  const stock = props.stock
  if (!stock) return []

  return [
    {
      label: t('inventory.stock_drawer.fields.created_at'),
      value: formatDateTime(
        stock.created_at,
        { dateStyle: 'medium', timeStyle: 'short' },
        t('inventory.stock_drawer.values.none'),
      ),
    },
    {
      label: t('inventory.stock_drawer.fields.updated_at'),
      value: formatDateTime(
        stock.updated_at,
        { dateStyle: 'medium', timeStyle: 'short' },
        t('inventory.stock_drawer.values.none'),
      ),
    },
    {
      label: t('inventory.stock_drawer.fields.serial_number'),
      value: toDisplayValue(props.materialDetail?.serial_number),
    },
    {
      label: t('inventory.stock_drawer.fields.order_number'),
      value: toDisplayValue(props.materialDetail?.order_number),
    },
    {
      label: t('inventory.stock_drawer.fields.lifetime_days'),
      value: toDisplayValue(props.materialDetail?.lifetime_days),
    },
  ]
})

const toggleMetadata = (): void => {
  isMetadataExpanded.value = !isMetadataExpanded.value
}
</script>

<template>
  <Transition name="inventory-drawer-fade">
    <div v-if="props.open" class="pointer-events-none fixed inset-0 z-[90]">
      <div class="pointer-events-auto absolute inset-x-0 top-24 bottom-0 bg-slate-900/35" @click="emit('close')" />

      <Transition name="inventory-drawer-slide">
        <aside
          v-if="props.stock"
          class="pointer-events-auto absolute top-24 right-0 bottom-0 flex w-full max-w-[46rem] flex-col border-l border-[var(--app-border)] bg-[var(--app-surface)] shadow-[0_16px_44px_rgba(15,23,42,0.18)]"
        >
          <header class="flex items-start justify-between gap-4 border-b border-[var(--app-border)] px-5 py-4">
            <div class="min-w-0 space-y-1">
              <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
                {{ t('inventory.stock_drawer.title') }}
              </p>
              <h2 class="truncate text-lg font-semibold text-slate-900">
                {{ props.stock.material.product_name }}
              </h2>
              <p class="text-xs text-slate-600">
                {{ t('inventory.stock_drawer.subtitle', { id: props.stock.id }) }}
              </p>
            </div>

            <UButton
              size="xs"
              variant="ghost"
              color="neutral"
              icon="i-heroicons-x-mark"
              :label="t('inventory.stock_drawer.actions.close')"
              @click="emit('close')"
            />
          </header>

          <div class="space-y-5 overflow-y-auto px-5 py-4">
            <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
              <div class="flex flex-wrap items-center gap-2">
                <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
                  {{ t('inventory.stock_drawer.sections.operational') }}
                </p>
                <UBadge :color="inventoryStatusColor" variant="soft">{{ inventoryStatusLabel }}</UBadge>
              </div>

              <div class="grid gap-2 sm:grid-cols-2">
                <div
                  v-for="field in operationalFields"
                  :key="field.label"
                  :class="[
                    'rounded-lg border border-slate-200 bg-slate-50 px-3 py-2',
                    field.wide ? 'sm:col-span-2' : '',
                  ]"
                >
                  <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">{{ field.label }}</p>
                  <p class="mt-1 text-sm text-slate-800">{{ field.value }}</p>
                </div>
              </div>
            </section>

            <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
              <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
                {{ t('inventory.stock_drawer.sections.unit_conversion') }}
              </p>

              <div class="grid gap-2 sm:grid-cols-3">
                <div
                  v-for="field in unitConversionFields"
                  :key="field.label"
                  class="rounded-lg border border-slate-200 bg-slate-50 px-3 py-2"
                >
                  <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">{{ field.label }}</p>
                  <p class="mt-1 text-sm text-slate-800">{{ field.value }}</p>
                </div>
              </div>
            </section>

            <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
              <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
                {{ t('inventory.stock_drawer.sections.material_identity') }}
              </p>

              <p v-if="props.isMaterialLoading" class="text-sm text-slate-500">
                {{ t('inventory.stock_workspace.loading') }}
              </p>

              <div v-else class="grid gap-2 sm:grid-cols-2">
                <div
                  v-for="field in materialIdentityFields"
                  :key="field.label"
                  :class="[
                    'rounded-lg border border-slate-200 bg-slate-50 px-3 py-2',
                    field.wide ? 'sm:col-span-2' : '',
                  ]"
                >
                  <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">{{ field.label }}</p>
                  <p class="mt-1 text-sm text-slate-800">{{ field.value }}</p>
                </div>
              </div>
            </section>

            <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
              <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
                {{ t('inventory.stock_drawer.sections.supplier') }}
              </p>

              <div class="grid gap-2 sm:grid-cols-2">
                <div
                  v-for="field in supplierFields"
                  :key="field.label"
                  :class="[
                    'rounded-lg border border-slate-200 bg-slate-50 px-3 py-2',
                    field.wide ? 'sm:col-span-2' : '',
                  ]"
                >
                  <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">{{ field.label }}</p>
                  <p class="mt-1 text-sm text-slate-800">{{ field.value }}</p>
                </div>
              </div>
            </section>

            <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
              <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
                {{ t('inventory.stock_drawer.sections.usage') }}
              </p>

              <div class="space-y-2 rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
                <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">
                  {{ t('inventory.stock_drawer.fields.recent_usages') }}
                </p>

                <p v-if="props.isUsagesLoading" class="text-sm text-slate-500">
                  {{ t('inventory.stock_workspace.loading') }}
                </p>
                <p v-else-if="sortedUsageEntries.length === 0" class="text-sm text-slate-600">
                  {{ t('inventory.stock_drawer.values.no_usage') }}
                </p>
                <ul v-else class="space-y-1.5 text-sm text-slate-700">
                  <li v-for="usage in sortedUsageEntries.slice(0, 8)" :key="usage.id">
                    <span class="font-medium">
                      {{
                        formatDateTime(
                          usage.used_at,
                          { dateStyle: 'medium', timeStyle: 'short' },
                          t('inventory.stock_drawer.values.none'),
                        )
                      }}
                    </span>
                    · {{ usage.quantity_used }} {{ usage.usage_unit.display_name }}
                    <span class="text-slate-500">
                      ({{ usage.project?.label || usage.project?.name || t('inventory.stock_drawer.values.none') }})
                    </span>
                  </li>
                </ul>
              </div>

              <div class="space-y-2 rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
                <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">
                  {{ t('inventory.stock_drawer.fields.linked_experiments_projects') }}
                </p>
                <p v-if="linkedExperimentProjectLabels.length === 0" class="text-sm text-slate-600">
                  {{ t('inventory.stock_drawer.values.no_linked_items') }}
                </p>
                <ul v-else class="list-disc space-y-1 pl-5 text-sm text-slate-700">
                  <li v-for="linkedLabel in linkedExperimentProjectLabels" :key="linkedLabel">
                    {{ linkedLabel }}
                  </li>
                </ul>
              </div>
            </section>

            <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
              <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
                {{ t('inventory.stock_drawer.sections.order') }}
              </p>

              <div class="space-y-2 rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
                <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">
                  {{ t('inventory.stock_drawer.fields.order_history') }}
                </p>

                <p v-if="props.isOrdersLoading" class="text-sm text-slate-500">
                  {{ t('inventory.stock_workspace.loading') }}
                </p>
                <p v-else-if="sortedOrderEntries.length === 0" class="text-sm text-slate-600">
                  {{ t('inventory.stock_drawer.values.no_order_history') }}
                </p>
                <ul v-else class="space-y-1.5 text-sm text-slate-700">
                  <li v-for="order in sortedOrderEntries.slice(0, 8)" :key="order.id">
                    <span class="font-medium">
                      {{
                        formatDateTime(
                          order.order_date,
                          { dateStyle: 'medium' },
                          t('inventory.stock_drawer.values.none'),
                        )
                      }}
                    </span>
                    · {{ order.amount }} {{ order.order_unit.display_name }}
                    <span class="text-slate-500"> ({{ order.status_label || order.status }}) </span>
                  </li>
                </ul>
              </div>

              <div class="space-y-2 rounded-lg border border-slate-200 bg-slate-50 px-3 py-2">
                <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">
                  {{ t('inventory.stock_drawer.fields.responsible_user') }}
                </p>
                <p v-if="responsibleUsers.length === 0" class="text-sm text-slate-600">
                  {{ t('inventory.stock_drawer.values.no_responsible_user') }}
                </p>
                <ul v-else class="list-disc space-y-1 pl-5 text-sm text-slate-700">
                  <li v-for="userLabel in responsibleUsers" :key="userLabel">
                    {{ userLabel }}
                  </li>
                </ul>
              </div>
            </section>

            <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
              <div class="flex items-center justify-between gap-3">
                <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
                  {{ t('inventory.stock_drawer.sections.metadata') }}
                </p>

                <UButton
                  size="xs"
                  variant="ghost"
                  color="neutral"
                  :label="
                    isMetadataExpanded
                      ? t('inventory.stock_drawer.metadata_toggle.hide')
                      : t('inventory.stock_drawer.metadata_toggle.show')
                  "
                  @click="toggleMetadata"
                />
              </div>

              <div v-if="isMetadataExpanded" class="grid gap-2 sm:grid-cols-2">
                <div
                  v-for="field in metadataFields"
                  :key="field.label"
                  class="rounded-lg border border-slate-200 bg-slate-50 px-3 py-2"
                >
                  <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">{{ field.label }}</p>
                  <p class="mt-1 text-sm text-slate-800">{{ field.value }}</p>
                </div>
              </div>
            </section>
          </div>
        </aside>
      </Transition>
    </div>
  </Transition>
</template>

<style scoped>
.inventory-drawer-fade-enter-active,
.inventory-drawer-fade-leave-active {
  transition: opacity 0.2s ease;
}

.inventory-drawer-fade-enter-from,
.inventory-drawer-fade-leave-to {
  opacity: 0;
}

.inventory-drawer-slide-enter-active,
.inventory-drawer-slide-leave-active {
  transition: transform 0.24s ease;
}

.inventory-drawer-slide-enter-from,
.inventory-drawer-slide-leave-to {
  transform: translateX(100%);
}
</style>
