<script setup lang="ts">
import { computed, onMounted, ref } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import WavesModalWrapper from '~/components/common/WavesModalWrapper.vue'
import { useManagementStore } from '~/stores/management'
import type { FileSystemItem } from '~/types/lab'
import type { NavigationTreeNode } from '~/types/navigation'
import { getErrorMessage } from '~/utils/errors'

const props = withDefaults(defineProps<{ filter?: string }>(), {
  filter: '',
})

const { t } = useI18n()
const toast = useToast()
const managementStore = useManagementStore()

const managementErrorMessage = ref<string | null>(null)
const isDeleteDialogOpen = ref(false)
const isUploadDialogOpen = ref(false)
const isFileContentDialogOpen = ref(false)
const uploadDirectoryPath = ref('')
const uploadFileValue = ref<File | null>(null)
const viewFilePath = ref('')
const viewFileContent = ref('')

onMounted(async () => {
  await initializeManagementTree()
})

const initializeManagementTree = async () => {
  managementErrorMessage.value = null

  try {
    await managementStore.initialize()
  } catch (err: unknown) {
    managementErrorMessage.value = getErrorMessage(err)
  }
}

/**
 * Converts one filesystem tree into UTree node format.
 *
 * Accepted data example:
 * - `{ type: 'directory', name: 'data', path: '/data', children: [{ type: 'file', name: 'a.txt', path: '/data/a.txt' }] }`
 *
 * Returned data example:
 * - `{ label: 'data', type: 'directory', path: '/data', children: [{ label: 'a.txt', type: 'file' }] }`
 */
const convertToTreeNode = (item: FileSystemItem): NavigationTreeNode => {
  const itemIcon = item.type === 'directory' ? 'i-heroicons-folder' : 'i-heroicons-document'
  const childNodes: NavigationTreeNode[] = []

  if (item.type === 'directory') {
    for (const child of item.children) {
      childNodes.push(convertToTreeNode(child))
    }
  }

  return {
    label: item.name,
    slot: 'management-path-node',
    icon: itemIcon,
    itemType: item.type,
    itemPath: item.path,
    draggable: true,
    onSelect: () => onPathNodeSelected(item.path),
    children: childNodes,
  }
}

const rootChildren = computed<NavigationTreeNode[]>(() => {
  const children: NavigationTreeNode[] = []
  const root = managementStore.dataDirectory

  if (root.type === 'directory') {
    for (const child of root.children) {
      children.push(convertToTreeNode(child))
    }
  }

  return children
})

const rawItems = computed<NavigationTreeNode[]>(() => [
  {
    label: t('navigation.management.root_label'),
    slot: 'management-root',
    icon: 'i-heroicons-cog-6-tooth',
    defaultExpanded: true,
    onSelect: () => navigateTo('/management'),
    children: rootChildren.value,
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

const getSlotItemLabel = (item: unknown): string => {
  if (typeof item !== 'object' || item === null || !('label' in item)) {
    return ''
  }

  const label = (item as { label?: unknown }).label
  return typeof label === 'string' ? label : ''
}

/**
 * Reads `itemPath` from one slot payload in a runtime-safe way.
 *
 * Accepted data examples:
 * - `{ itemPath: '/data' }`
 * - `{ itemPath: '/data/file.txt' }`
 *
 * Returned data examples:
 * - `'/data'`
 * - `''` when path is missing
 */
const getSlotItemPath = (item: unknown): string => {
  if (typeof item !== 'object' || item === null || !('itemPath' in item)) {
    return ''
  }

  const path = (item as { itemPath?: unknown }).itemPath
  return typeof path === 'string' ? path : ''
}

const onPathNodeSelected = (path: string) => {
  managementStore.selectedPath = path
  managementStore.selectedPaths.push(path)
}

const openDeleteDialog = (path: string) => {
  managementStore.selectedPath = path
  isDeleteDialogOpen.value = true
}

const confirmDelete = async () => {
  const path = managementStore.selectedPath
  if (path.trim() === '') {
    return
  }

  try {
    await managementStore.deleteFile(path)
    isDeleteDialogOpen.value = false

    toast.add({
      title: t('navigation.management.toasts.delete_success'),
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: t('navigation.management.toasts.delete_failed'),
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  }
}

const getFileNameFromPath = (path: string): string => {
  const parts = path.split('/')
  const fileName = parts[parts.length - 1]
  return fileName && fileName !== '' ? fileName : path
}

const downloadFile = async (path: string) => {
  try {
    const fileBlob = await managementStore.downloadFile(path)

    const downloadUrl = window.URL.createObjectURL(fileBlob)
    const link = document.createElement('a')
    link.href = downloadUrl
    link.setAttribute('download', getFileNameFromPath(path))
    document.body.appendChild(link)
    link.click()
    document.body.removeChild(link)
    window.URL.revokeObjectURL(downloadUrl)
  } catch (err: unknown) {
    toast.add({
      title: t('navigation.management.toasts.download_failed'),
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  }
}

const openUploadDialog = (path: string) => {
  uploadDirectoryPath.value = path
  uploadFileValue.value = null
  isUploadDialogOpen.value = true
}

const onUploadFileInput = (event: Event) => {
  const target = event.target as HTMLInputElement | null
  if (!target || !target.files || target.files.length === 0) {
    uploadFileValue.value = null
    return
  }

  const selectedFile = target.files[0]
  uploadFileValue.value = selectedFile ?? null
}

const confirmUpload = async () => {
  const selectedFile = uploadFileValue.value
  if (!selectedFile) {
    toast.add({
      title: t('navigation.management.toasts.upload_file_required'),
      color: 'warning',
      duration: 2500,
    })
    return
  }

  try {
    await managementStore.uploadFile(uploadDirectoryPath.value, selectedFile)
    isUploadDialogOpen.value = false

    toast.add({
      title: t('navigation.management.toasts.upload_success'),
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: t('navigation.management.toasts.upload_failed'),
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  }
}

const openFileContentDialog = async (path: string) => {
  try {
    const content = await managementStore.getFileContent(path)
    viewFilePath.value = path
    viewFileContent.value = content
    isFileContentDialogOpen.value = true
  } catch (err: unknown) {
    toast.add({
      title: t('navigation.management.toasts.file_content_failed'),
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  }
}

const handleDragStart = (path: string, event: DragEvent) => {
  event.dataTransfer?.setData('text/plain', path)
}

const canShowContentAction = (label: string): boolean => {
  return !/\.(png|jpe?g)$/i.test(label)
}

const getPathContextMenuItems = (path: string, itemType: string, label: string) => {
  const firstGroup = [
    {
      label: t('navigation.management.actions.delete_item'),
      icon: 'i-heroicons-trash',
      onSelect: () => openDeleteDialog(path),
    },
  ]

  const secondGroup: Array<{ label: string; icon: string; onSelect: () => void }> = []

  if (itemType === 'file') {
    secondGroup.push({
      label: t('navigation.management.actions.download_file'),
      icon: 'i-heroicons-arrow-down-tray',
      onSelect: () => downloadFile(path),
    })

    if (canShowContentAction(label)) {
      secondGroup.push({
        label: t('navigation.management.actions.view_file_content'),
        icon: 'i-heroicons-document-magnifying-glass',
        onSelect: () => {
          void openFileContentDialog(path)
        },
      })
    }
  }

  if (itemType === 'directory') {
    secondGroup.push({
      label: t('navigation.management.actions.upload_file'),
      icon: 'i-heroicons-arrow-up-tray',
      onSelect: () => openUploadDialog(path),
    })
  }

  const menuItems: Array<Array<{ label: string; icon: string; onSelect: () => void }>> = []
  menuItems.push(firstGroup)
  if (secondGroup.length > 0) {
    menuItems.push(secondGroup)
  }

  return menuItems
}

const getPathContextMenuItemsFromSlotItem = (item: unknown) => {
  if (typeof item !== 'object' || item === null) {
    return []
  }

  if (!('itemPath' in item) || !('itemType' in item) || !('label' in item)) {
    return []
  }

  const path = (item as { itemPath?: unknown }).itemPath
  const itemType = (item as { itemType?: unknown }).itemType
  const label = (item as { label?: unknown }).label

  if (typeof path !== 'string' || typeof itemType !== 'string' || typeof label !== 'string') {
    return []
  }

  return getPathContextMenuItems(path, itemType, label)
}
</script>

<template>
  <section class="space-y-1">
    <p class="px-2 text-xs font-semibold tracking-wide text-slate-500 uppercase">
      {{ t('navigation.management.title') }}
    </p>
    <p v-if="managementStore.isLoadingDirectoryContent" class="px-2 text-xs text-slate-500">
      {{ t('navigation.management.loading') }}
    </p>
    <p v-else-if="managementErrorMessage" class="px-2 text-xs text-red-600">
      {{ managementErrorMessage }}
    </p>
    <p v-else-if="rootChildren.length === 0" class="px-2 text-xs text-slate-500">
      {{ t('navigation.management.empty') }}
    </p>

    <UTree v-else :items="items" :ui="{ link: 'cursor-pointer hover:text-primary transition-colors' }">
      <template #management-root-label="{ item }">
        <span class="inline-block w-full">{{ getSlotItemLabel(item) }}</span>
      </template>

      <template #management-path-node-label="{ item }">
        <UContextMenu :items="getPathContextMenuItemsFromSlotItem(item)">
          <span
            class="fileSystemItem inline-block w-full"
            draggable="true"
            @dragstart="handleDragStart(getSlotItemPath(item), $event)"
          >
            {{ getSlotItemLabel(item) }}
          </span>
        </UContextMenu>
      </template>
    </UTree>

    <WavesModalWrapper
      :open="isDeleteDialogOpen"
      :title="t('navigation.management.delete_modal.title')"
      :description="t('navigation.management.delete_modal.description', { path: managementStore.selectedPath })"
      @update:open="isDeleteDialogOpen = $event"
    >
      <template #footer>
        <BaseButton
          :label="t('common.actions.cancel')"
          :on-click="() => (isDeleteDialogOpen = false)"
          variant="secondary"
          size="sm"
          width="auto"
          :disabled="managementStore.isDeletingFile"
        />
        <BaseButton
          :label="t('common.actions.delete')"
          :on-click="confirmDelete"
          variant="danger"
          size="sm"
          width="auto"
          :loading="managementStore.isDeletingFile"
        />
      </template>
    </WavesModalWrapper>

    <WavesModalWrapper
      :open="isUploadDialogOpen"
      :title="t('navigation.management.upload_modal.title')"
      :description="t('navigation.management.upload_modal.description', { path: uploadDirectoryPath })"
      @update:open="isUploadDialogOpen = $event"
    >
      <template #body>
        <input
          type="file"
          class="w-full rounded-xl border border-black/15 bg-white/70 px-4 py-3 text-sm text-slate-700"
          @change="onUploadFileInput"
        />
      </template>

      <template #footer>
        <BaseButton
          :label="t('common.actions.cancel')"
          :on-click="() => (isUploadDialogOpen = false)"
          variant="secondary"
          size="sm"
          width="auto"
          :disabled="managementStore.isUploadingFile"
        />
        <BaseButton
          :label="t('navigation.management.upload_modal.upload_button')"
          :on-click="confirmUpload"
          variant="primary"
          size="sm"
          width="auto"
          :loading="managementStore.isUploadingFile"
        />
      </template>
    </WavesModalWrapper>

    <UModal
      :open="isFileContentDialogOpen"
      :dismissible="true"
      fullscreen
      @update:open="isFileContentDialogOpen = $event"
    >
      <template #body>
        <div class="h-full bg-slate-900 text-white">
          <div class="flex items-center justify-between border-b border-slate-700 px-4 py-3">
            <p class="truncate text-sm font-semibold">{{ viewFilePath }}</p>
            <UButton
              size="xs"
              variant="outline"
              color="neutral"
              icon="i-heroicons-x-mark"
              :label="t('common.actions.close')"
              @click="isFileContentDialogOpen = false"
            />
          </div>

          <div class="h-[calc(100%-52px)] overflow-auto p-4">
            <pre class="text-xs leading-5 whitespace-pre-wrap">{{ viewFileContent }}</pre>
          </div>
        </div>
      </template>
    </UModal>
  </section>
</template>

<style scoped>
.fileSystemItem:hover {
  color: #1d4ed8;
}
</style>
