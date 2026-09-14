<script setup lang="ts">
import { computed, ref } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import { useCompoundLibraryStore } from '~/stores/compoundLibraries'
import { usePlateStore } from '~/stores/plates'
import { getErrorMessage } from '~/utils/errors'

const { t } = useI18n()
const toast = useToast()
const plateStore = usePlateStore()
const compoundLibraryStore = useCompoundLibraryStore()

const isConfirmOpen = ref(false)

// An unset value (null in the database) counts as "not archived", as in the navigation tree.
const isArchived = computed(() => Boolean(plateStore.currentPlate?.archived))
const barcode = computed(() => plateStore.currentPlate?.barcode ?? '')

const openConfirmation = () => {
  if (plateStore.isArchivingPlate) return
  isConfirmOpen.value = true
}

const closeConfirmation = () => {
  if (plateStore.isArchivingPlate) return
  isConfirmOpen.value = false
}

const onConfirmOpenChange = (isOpen: boolean) => {
  if (isOpen) {
    openConfirmation()
  } else {
    closeConfirmation()
  }
}

/**
 * Switches the plate between archived and active after the user confirmed it.
 */
const confirmChange = async () => {
  const shouldArchive = !isArchived.value

  try {
    const result = await plateStore.setCurrentPlateArchived(shouldArchive)
    // The navigation tree keeps its own copy of the library plates.
    compoundLibraryStore.setPlateArchived(result.id, result.archived)
    isConfirmOpen.value = false
    toast.add({
      title: shouldArchive ? t('plates.page.archive.archived_toast') : t('plates.page.archive.unarchived_toast'),
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: t('plates.page.archive.failed_toast'),
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  }
}
</script>

<template>
  <div class="flex shrink-0 items-center gap-3">
    <span v-if="isArchived" class="rounded-full bg-slate-200 px-3 py-1 text-xs font-semibold text-slate-700">
      {{ t('plates.page.archive.archived_badge') }}
    </span>

    <BaseButton
      :label="isArchived ? t('plates.page.archive.unarchive_button') : t('plates.page.archive.archive_button')"
      :on-click="openConfirmation"
      variant="secondary"
      size="sm"
      width="auto"
      :disabled="plateStore.isArchivingPlate"
    />

    <!-- A plain confirmation card, the same style as archiving a stock item in the inventory. -->
    <UModal
      :open="isConfirmOpen"
      :title="isArchived ? t('plates.page.archive.unarchive_title') : t('plates.page.archive.archive_title')"
      :description="
        isArchived
          ? t('plates.page.archive.unarchive_description', { barcode })
          : t('plates.page.archive.archive_description', { barcode })
      "
      :dismissible="!plateStore.isArchivingPlate"
      class="w-full sm:max-w-lg"
      :ui="{ content: 'rounded-2xl bg-white shadow-md' }"
      @update:open="onConfirmOpenChange"
    >
      <template #footer>
        <div class="flex w-full justify-end gap-2">
          <UButton
            variant="ghost"
            color="neutral"
            :label="t('common.actions.cancel')"
            :disabled="plateStore.isArchivingPlate"
            @click="closeConfirmation"
          />
          <UButton
            color="warning"
            :label="isArchived ? t('plates.page.archive.unarchive_confirm') : t('plates.page.archive.archive_confirm')"
            :loading="plateStore.isArchivingPlate"
            :disabled="plateStore.isArchivingPlate"
            @click="confirmChange"
          />
        </div>
      </template>
    </UModal>
  </div>
</template>
