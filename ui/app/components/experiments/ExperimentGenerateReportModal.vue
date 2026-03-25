<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import WavesModalWrapper from '~/components/common/WavesModalWrapper.vue'
import { useReportStore } from '~/stores/reports'
import { getErrorMessage } from '~/utils/errors'

const props = defineProps<{
  open: boolean
  experimentName: string
  labels: string[]
  selectedPos: string | null
  selectedNeg: string | null
}>()

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
  (e: 'generated'): void
}>()

const { t } = useI18n()
const toast = useToast()
const reportStore = useReportStore()

const selectedLabel = ref('')
const selectedNotebookPath = ref('')

watch(
  () => props.open,
  async (isOpen) => {
    if (!isOpen) {
      return
    }

    selectedLabel.value = ''
    selectedNotebookPath.value = ''

    if (props.labels.length > 0) {
      const firstLabel = props.labels[0]
      if (firstLabel) {
        selectedLabel.value = firstLabel
      }
    }

    await reportStore.fetchInputNotebookOptions()
  },
  { immediate: true },
)

/**
 * Checks whether both controls were saved in the "Show results" section.
 *
 * Returned data examples:
 * - `true` when `selectedPos = 'P'` and `selectedNeg = 'N'`
 * - `false` when one control is missing
 */
const hasSavedControls = computed(() => {
  if (!props.selectedPos) {
    return false
  }

  if (!props.selectedNeg) {
    return false
  }

  return true
})

/**
 * Enables the generate button only when all required fields are set.
 *
 * Required values:
 * - positive and negative control
 * - notebook template path
 * - measurement label
 */
const canGenerate = computed(() => {
  if (!hasSavedControls.value) {
    return false
  }

  if (selectedNotebookPath.value.trim() === '') {
    return false
  }

  if (selectedLabel.value.trim() === '') {
    return false
  }

  if (reportStore.isGeneratingReport) {
    return false
  }

  return true
})

const close = () => {
  emit('update:open', false)
}

/**
 * Converts one notebook path into a short option label.
 *
 * Accepted input example:
 * - `'/notebooks/input/default_report.ipynb'`
 *
 * Returned data example:
 * - `'default_report.ipynb'`
 */
const formatNotebookOptionLabel = (notebookPath: string): string => {
  return reportStore.getFileNameFromPath(notebookPath)
}

const submit = async () => {
  if (!hasSavedControls.value) {
    toast.add({
      title: t('experiments.reports.controls_required'),
      color: 'error',
      duration: 3500,
    })
    return
  }

  if (selectedNotebookPath.value.trim() === '') {
    toast.add({
      title: t('experiments.reports.select_template_error'),
      color: 'error',
      duration: 3500,
    })
    return
  }

  if (selectedLabel.value.trim() === '') {
    toast.add({
      title: t('experiments.reports.select_label_error'),
      color: 'error',
      duration: 3500,
    })
    return
  }

  if (!props.selectedPos || !props.selectedNeg) {
    return
  }

  try {
    await reportStore.generateReport({
      experiment: props.experimentName,
      label: selectedLabel.value,
      notebook_path: selectedNotebookPath.value,
      selected_pos: props.selectedPos,
      selected_neg: props.selectedNeg,
    })

    emit('generated')
    close()

    toast.add({
      title: t('experiments.reports.generated_success'),
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: t('experiments.reports.generate_failed'),
      description: getErrorMessage(err),
      color: 'error',
      duration: 4500,
    })
  }
}
</script>

<template>
  <WavesModalWrapper
    :open="props.open"
    :title="t('experiments.reports.modal_title')"
    :description="t('experiments.reports.modal_description')"
    :dismissible="!reportStore.isGeneratingReport"
    modal-class="w-full sm:max-w-lg"
    body-container-class="w-full px-8 pt-10"
    @update:open="emit('update:open', $event)"
  >
    <template #body>
      <div class="space-y-4">
        <div>
          <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
            {{ t('experiments.reports.notebook_template_label') }}
          </label>

          <select
            v-model="selectedNotebookPath"
            class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-blue-300"
          >
            <option value="">
              {{ t('experiments.reports.notebook_template_placeholder') }}
            </option>
            <option
              v-for="notebookPath in reportStore.inputNotebookOptions"
              :key="`report-template-${notebookPath}`"
              :value="notebookPath"
            >
              {{ formatNotebookOptionLabel(notebookPath) }}
            </option>
          </select>
        </div>

        <div>
          <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
            {{ t('experiments.reports.measurement_label') }}
          </label>

          <select
            v-model="selectedLabel"
            class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-blue-300"
          >
            <option value="">
              {{ t('experiments.reports.measurement_placeholder') }}
            </option>
            <option v-for="label in props.labels" :key="`report-label-${label}`" :value="label">
              {{ label }}
            </option>
          </select>
        </div>

        <p v-if="!hasSavedControls" class="text-sm text-red-600">
          {{ t('experiments.reports.controls_required') }}
        </p>
      </div>
    </template>

    <template #footer>
      <BaseButton
        :label="t('common.actions.cancel')"
        :on-click="close"
        variant="secondary"
        size="sm"
        width="auto"
        :disabled="reportStore.isGeneratingReport"
      />
      <BaseButton
        :label="t('experiments.reports.generate_button')"
        :on-click="submit"
        variant="primary"
        size="sm"
        width="auto"
        :loading="reportStore.isGeneratingReport"
        :disabled="!canGenerate"
      />
    </template>
  </WavesModalWrapper>
</template>
