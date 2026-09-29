<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import WavesModalWrapper from '~/components/common/WavesModalWrapper.vue'
import { useAnalysisStore } from '~/stores/analysis'
import type { AnalysisType } from '~/types/analysis'
import type { CommandMessage } from '~/types/management'

const props = defineProps<{
  open: boolean
  experimentId: number
  labels: string[]
}>()

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
}>()

const { t } = useI18n()
const analysisStore = useAnalysisStore()

const analysisType = ref<AnalysisType>('single')
const selectedLabel = ref('')
const conditionYes = ref('')
const conditionNo = ref('')

watch(
  () => props.open,
  (isOpen) => {
    // The page is kept when another experiment is opened, so its labels may differ
    if (isOpen && !props.labels.includes(selectedLabel.value)) {
      selectedLabel.value = props.labels[0] ?? ''
    }
  },
  { immediate: true },
)

// The messages of the last run belong to its experiment only
const isThisExperiment = computed(() => analysisStore.experimentId === props.experimentId)
const isRunningElsewhere = computed(() => analysisStore.isRunning && !isThisExperiment.value)

const canStart = computed(() => {
  if (analysisStore.isRunning) return false
  if (selectedLabel.value === '') return false
  if (analysisType.value === 'selectivity') {
    return conditionYes.value.trim() !== '' && conditionNo.value.trim() !== ''
  }
  return true
})

const start = async () => {
  const settings =
    analysisType.value === 'selectivity'
      ? { condi_yes: conditionYes.value.trim(), condi_no: conditionNo.value.trim() }
      : {}
  await analysisStore.startAnalysis({
    experiment_id: props.experimentId,
    label: selectedLabel.value,
    analysis_type: analysisType.value,
    settings,
  })
}

const close = () => {
  emit('update:open', false)
}

// The summary above the messages, once the analysis has ended
const summary = computed(() => {
  if (!isThisExperiment.value) {
    return null
  }
  if (analysisStore.status === 'failed') {
    return { color: 'error' as const, icon: 'i-lucide-circle-x', title: t('experiments.analysis.failed') }
  }
  if (analysisStore.status === 'completed') {
    return { color: 'success' as const, icon: 'i-lucide-circle-check', title: t('experiments.analysis.completed') }
  }
  return null
})

const messageClass = (message: CommandMessage): string => {
  if (message.level === 'error') return 'bg-red-50 text-red-800'
  if (message.level === 'warning') return 'bg-amber-50 text-amber-800'
  if (message.level === 'success') return 'text-green-700'
  return 'text-slate-700'
}

const fieldClass =
  'w-full rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-blue-300'
</script>

<template>
  <WavesModalWrapper
    :open="props.open"
    :title="t('experiments.analysis.modal_title')"
    :description="t('experiments.analysis.modal_description')"
    modal-class="w-full sm:max-w-2xl"
    body-container-class="w-full px-8 pt-10"
    @update:open="emit('update:open', $event)"
  >
    <template #body>
      <div class="space-y-4">
        <div>
          <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
            {{ t('experiments.analysis.type') }}
          </label>
          <select v-model="analysisType" :class="[fieldClass, 'cursor-pointer']" :disabled="analysisStore.isRunning">
            <option value="single">{{ t('experiments.analysis.type_single') }}</option>
            <option value="selectivity">{{ t('experiments.analysis.type_selectivity') }}</option>
          </select>
        </div>

        <div>
          <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
            {{ t('experiments.analysis.measurement_label') }}
          </label>
          <select v-model="selectedLabel" :class="[fieldClass, 'cursor-pointer']" :disabled="analysisStore.isRunning">
            <option v-for="label in props.labels" :key="`analysis-label-${label}`" :value="label">
              {{ label }}
            </option>
          </select>
        </div>

        <div v-if="analysisType === 'selectivity'" class="grid grid-cols-2 gap-3">
          <div>
            <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
              {{ t('experiments.analysis.condition_yes') }}
            </label>
            <input v-model="conditionYes" type="text" :class="fieldClass" :disabled="analysisStore.isRunning" />
          </div>
          <div>
            <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
              {{ t('experiments.analysis.condition_no') }}
            </label>
            <input v-model="conditionNo" type="text" :class="fieldClass" :disabled="analysisStore.isRunning" />
          </div>
        </div>

        <p v-if="isRunningElsewhere" class="text-sm text-amber-700">
          {{ t('experiments.analysis.running_elsewhere', { id: analysisStore.experimentId }) }}
        </p>

        <div
          v-if="isThisExperiment && (analysisStore.messages.length > 0 || analysisStore.isRunning)"
          class="space-y-2"
        >
          <UAlert v-if="summary" :color="summary.color" :icon="summary.icon" :title="summary.title" variant="subtle" />
          <div class="max-h-80 overflow-auto rounded-xl border border-slate-200 bg-white p-2">
            <p
              v-for="(message, index) in analysisStore.messages"
              :key="index"
              class="rounded px-2 py-0.5 font-mono text-xs whitespace-pre-wrap"
              :class="messageClass(message)"
              v-text="message.text"
            />
            <p v-if="analysisStore.isRunning" class="px-2 py-0.5 text-xs text-slate-500">
              {{ t('experiments.analysis.running') }}
            </p>
          </div>
        </div>
      </div>
    </template>

    <template #footer>
      <BaseButton :label="t('common.actions.close')" :on-click="close" variant="secondary" size="sm" width="auto" />
      <BaseButton
        :label="t('experiments.analysis.start_button')"
        :on-click="start"
        variant="primary"
        size="sm"
        width="auto"
        :loading="analysisStore.isRunning"
        :disabled="!canStart"
      />
    </template>
  </WavesModalWrapper>
</template>
