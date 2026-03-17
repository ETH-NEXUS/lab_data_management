<script setup lang="ts">
import { computed, onMounted, ref } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import BaseField from '~/components/common/BaseField.vue'
import WavesModalWrapper from '~/components/common/WavesModalWrapper.vue'
import { useCompoundLibraryStore } from '~/stores/compoundLibraries'
import type { CompoundLibrary, Plate } from '~/types/lab'
import type { NavigationTreeNode } from '~/types/navigation'
import { getErrorMessage } from '~/utils/errors'

const props = withDefaults(defineProps<{ filter?: string }>(), {
  filter: '',
})

const { t } = useI18n()
const toast = useToast()
const compoundLibraryStore = useCompoundLibraryStore()

const librariesErrorMessage = ref<string | null>(null)
const isCreatePlateModalOpen = ref(false)
const selectedLibraryId = ref<number | null>(null)
const selectedLibraryName = ref('')
const newPlateBarcode = ref('')

onMounted(async () => {
  await loadLibraries()
})

/**
 * Loads all compound libraries used in this navigation branch.
 *
 * Loaded data example:
 * - `[{ id: 3, name: 'Kinase Library', file_name: 'kinase.csv', plates: [{ id: 44, barcode: 'LIB-001' }] }]`
 */
const loadLibraries = async () => {
  librariesErrorMessage.value = null

  try {
    await compoundLibraryStore.fetchLibraries()
  } catch (err: unknown) {
    librariesErrorMessage.value = getErrorMessage(err)
  }
}

/**
 * Formats plate dimension from different backend shapes into one readable label.
 *
 * Accepted data examples:
 * - `dimension = '384'`
 * - `dimension = { id: 1, name: '384', cols: 24, rows: 16 }`
 *
 * Returned data examples:
 * - `'384'`
 * - `'No dimension'`
 */
const formatPlateDimensionLabel = (dimension: unknown): string => {
  if (typeof dimension === 'string' || typeof dimension === 'number') {
    return String(dimension)
  }

  if (typeof dimension === 'object' && dimension !== null && 'name' in dimension) {
    const name = (dimension as { name?: unknown }).name
    if (typeof name === 'string' && name.trim() !== '') {
      return name
    }
  }

  return t('navigation.libraries.no_dimension')
}

/**
 * Creates sorted visible plate nodes for one library.
 *
 * Rules copied from old UI:
 * - archived plates are hidden
 * - plates are sorted by barcode
 * - `empty_wells` status shows warning symbol
 */
const mapLibraryPlateNodes = (library: CompoundLibrary): NavigationTreeNode[] => {
  const visiblePlates: Plate[] = []
  const libraryPlates = library.plates ?? []

  for (const plate of libraryPlates) {
    if (plate.archived) {
      continue
    }
    visiblePlates.push(plate)
  }

  visiblePlates.sort((leftPlate, rightPlate) => leftPlate.barcode.localeCompare(rightPlate.barcode))

  const nodes: NavigationTreeNode[] = []
  for (const plate of visiblePlates) {
    const plateLabel = plate.barcode || t('navigation.libraries.plate_fallback', { id: plate.id })
    const dimensionLabel = formatPlateDimensionLabel(plate.dimension)
    const hasWarning = plate.status === 'empty_wells'
    const warningLabel = hasWarning ? ' ⚠️' : ''
    const routePlateId = plate.barcode || String(plate.id)

    nodes.push({
      label: `${plateLabel} (${dimensionLabel})${warningLabel}`,
      icon: 'i-heroicons-squares-2x2',
      onSelect: () => navigateTo(`/plates/${encodeURIComponent(routePlateId)}`),
    })
  }

  return nodes
}

/**
 * Creates first-level compound-library nodes for UTree.
 *
 * Returned data example:
 * - `[{ label: 'Kinase Library', libraryId: 3, children: [{ label: 'LIB-001 (384)' }] }]`
 */
const libraryNodes = computed<NavigationTreeNode[]>(() => {
  const nodes: NavigationTreeNode[] = []

  for (const library of compoundLibraryStore.libraries) {
    const libraryLabel = library.name || t('navigation.libraries.library_fallback', { id: library.id })

    nodes.push({
      label: libraryLabel,
      slot: 'compound-library-node',
      libraryId: library.id,
      libraryName: libraryLabel,
      icon: 'i-heroicons-beaker',
      children: mapLibraryPlateNodes(library),
    })
  }

  return nodes
})

const rawItems = computed<NavigationTreeNode[]>(() => [
  {
    label: t('navigation.libraries.root_label'),
    slot: 'libraries-root',
    icon: 'i-heroicons-rectangle-stack',
    defaultExpanded: true,
    children: libraryNodes.value,
  },
])

const filterTree = (items: NavigationTreeNode[], query: string): NavigationTreeNode[] => {
  const q = query.trim().toLowerCase()
  if (!q) {
    return items
  }

  const filteredItems: NavigationTreeNode[] = []

  for (const item of items) {
    const filteredChildren = item.children ? filterTree(item.children, q) : []
    const selfMatches = item.label.toLowerCase().includes(q)
    const hasMatchingChildren = filteredChildren.length > 0

    if (!selfMatches && !hasMatchingChildren) {
      continue
    }

    filteredItems.push({
      ...item,
      children: filteredChildren,
      defaultExpanded: true,
    })
  }

  return filteredItems
}

const items = computed<NavigationTreeNode[]>(() => filterTree(rawItems.value, props.filter))

/**
 * Runtime-safe extraction of label text from UTree slot payload.
 *
 * Accepted data examples:
 * - `{ label: 'Libraries' }`
 * - `{ label: 'Kinase Library', libraryId: 3 }`
 */
const getSlotItemLabel = (item: unknown): string => {
  if (typeof item !== 'object' || item === null || !('label' in item)) {
    return ''
  }

  const label = (item as { label?: unknown }).label
  return typeof label === 'string' ? label : ''
}

const openCreatePlateModal = (libraryId: number, libraryName: string) => {
  selectedLibraryId.value = libraryId
  selectedLibraryName.value = libraryName
  newPlateBarcode.value = ''
  isCreatePlateModalOpen.value = true
}

const closeCreatePlateModal = () => {
  isCreatePlateModalOpen.value = false
  selectedLibraryId.value = null
  selectedLibraryName.value = ''
  newPlateBarcode.value = ''
}

const getLibraryContextMenuItems = (libraryId: number, libraryName: string) => [
  [
    {
      label: t('navigation.libraries.actions.new_plate'),
      icon: 'i-heroicons-plus-circle',
      onSelect: () => openCreatePlateModal(libraryId, libraryName),
    },
  ],
]

/**
 * Converts one UTree slot payload into context-menu actions.
 *
 * Accepted data example:
 * - `{ label: 'Kinase Library', libraryId: 3, libraryName: 'Kinase Library' }`
 */
const getLibraryContextMenuItemsFromSlotItem = (item: unknown) => {
  if (typeof item !== 'object' || item === null) {
    return []
  }

  if (!('libraryId' in item) || !('libraryName' in item)) {
    return []
  }

  const libraryId = (item as { libraryId?: unknown }).libraryId
  const libraryName = (item as { libraryName?: unknown }).libraryName

  if (typeof libraryId !== 'number' || typeof libraryName !== 'string') {
    return []
  }

  return getLibraryContextMenuItems(libraryId, libraryName)
}

const createPlateForLibrary = async () => {
  if (!selectedLibraryId.value) {
    return
  }

  const trimmedBarcode = newPlateBarcode.value.trim()
  if (trimmedBarcode === '') {
    toast.add({
      title: t('navigation.libraries.toasts.barcode_required'),
      color: 'warning',
      duration: 2500,
    })
    return
  }

  try {
    await compoundLibraryStore.addPlateToLibrary(selectedLibraryId.value, trimmedBarcode)
    closeCreatePlateModal()

    toast.add({
      title: t('navigation.libraries.toasts.created'),
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: t('navigation.libraries.toasts.create_failed'),
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  }
}
</script>

<template>
  <section class="space-y-1">
    <p class="px-2 text-xs font-semibold tracking-wide text-slate-500 uppercase">
      {{ t('navigation.libraries.title') }}
    </p>
    <p v-if="compoundLibraryStore.isLoadingLibraries" class="px-2 text-xs text-slate-500">
      {{ t('navigation.libraries.loading') }}
    </p>
    <p v-else-if="librariesErrorMessage" class="px-2 text-xs text-red-600">
      {{ librariesErrorMessage }}
    </p>
    <p v-else-if="libraryNodes.length === 0" class="px-2 text-xs text-slate-500">
      {{ t('navigation.libraries.empty') }}
    </p>

    <UTree v-else :items="items" :ui="{ link: 'cursor-pointer hover:text-primary transition-colors' }">
      <template #libraries-root-label="{ item }">
        <span class="inline-block w-full">{{ getSlotItemLabel(item) }}</span>
      </template>

      <template #compound-library-node-label="{ item }">
        <UContextMenu :items="getLibraryContextMenuItemsFromSlotItem(item)">
          <span class="inline-block w-full">{{ getSlotItemLabel(item) }}</span>
        </UContextMenu>
      </template>
    </UTree>

    <WavesModalWrapper
      :open="isCreatePlateModalOpen"
      :title="t('navigation.libraries.modal.title')"
      :description="t('navigation.libraries.modal.description', { libraryName: selectedLibraryName })"
      @update:open="isCreatePlateModalOpen = $event"
    >
      <template #body>
        <BaseField
          v-model="newPlateBarcode"
          :label="t('navigation.libraries.modal.barcode_label')"
          :placeholder="t('navigation.libraries.modal.barcode_placeholder')"
          input-class="w-full rounded-full bg-white/70 px-4 py-3 shadow outline-none focus:ring-2 focus:ring-lime-500"
        />
      </template>

      <template #footer>
        <BaseButton
          :label="t('common.actions.cancel')"
          :on-click="closeCreatePlateModal"
          variant="secondary"
          size="sm"
          width="auto"
          :disabled="compoundLibraryStore.isAddingLibraryPlate"
        />
        <BaseButton
          :label="t('navigation.libraries.modal.create_button')"
          :on-click="createPlateForLibrary"
          variant="primary"
          size="sm"
          width="auto"
          :loading="compoundLibraryStore.isAddingLibraryPlate"
        />
      </template>
    </WavesModalWrapper>
  </section>
</template>
