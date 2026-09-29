<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import CommandOutputLog from '~/components/common/CommandOutputLog.vue'
import WavesModalWrapper from '~/components/common/WavesModalWrapper.vue'
import { useAnalysisStore } from '~/stores/analysis'
import { usePlateViewStore } from '~/stores/plateView'
import type { AnalysisType } from '~/types/analysis'
import type { ExperimentDetails } from '~/types/lab'

const props = defineProps<{
  open: boolean
  experimentId: number
  labels: string[]
  // The well types of each label are its keys, e.g. { Lum: { C: {...}, N1: {...}, P1: {...} } }
  stats: ExperimentDetails['stats']
}>()

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
}>()

const { t } = useI18n()
const analysisStore = useAnalysisStore()
const plateViewStore = usePlateViewStore()

const analysisType = ref<AnalysisType>('single')
const selectedLabel = ref('')
const conditionYes = ref('')
const conditionNo = ref('')
const positiveControl = ref('')
const negativeControl = ref('')

// The well types of the chosen measurement, e.g. ['C', 'N1', 'P1']
const wellTypes = computed(() => Object.keys(props.stats[selectedLabel.value] ?? {}).sort())

/**
 * The control that is chosen first: the one saved under "Show results", else the
 * type called exactly 'P' (or 'N'), else the first numbered one, e.g. 'P1'.
 *
 * Returned data examples: 'P', 'P1', '' (no such well type)
 */
const defaultControl = (saved: string | null | undefined, letter: 'P' | 'N'): string => {
  if (saved && wellTypes.value.includes(saved)) return saved
  if (wellTypes.value.includes(letter)) return letter
  return wellTypes.value.find((wellType) => wellType.startsWith(letter)) ?? ''
}

const chooseDefaultControls = () => {
  const saved = plateViewStore.getExperimentControls(props.experimentId)
  if (!wellTypes.value.includes(positiveControl.value)) {
    positiveControl.value = defaultControl(saved?.pos, 'P')
  }
  if (!wellTypes.value.includes(negativeControl.value)) {
    negativeControl.value = defaultControl(saved?.neg, 'N')
  }
}

watch(
  () => props.open,
  (isOpen) => {
    // The page is kept when another experiment is opened, so its labels may differ
    if (isOpen && !props.labels.includes(selectedLabel.value)) {
      selectedLabel.value = props.labels[0] ?? ''
    }
    if (isOpen) chooseDefaultControls()
  },
  { immediate: true },
)

watch(selectedLabel, chooseDefaultControls)

// The page is kept when another experiment is opened: its controls are chosen anew
watch(
  () => props.experimentId,
  () => {
    positiveControl.value = ''
    negativeControl.value = ''
  },
)

// The messages of the last run belong to its experiment only
const isThisExperiment = computed(() => analysisStore.experimentId === props.experimentId)
const isRunningElsewhere = computed(() => analysisStore.isRunning && !isThisExperiment.value)

const canStart = computed(() => {
  if (analysisStore.isRunning) return false
  if (selectedLabel.value === '') return false
  if (positiveControl.value === '' || negativeControl.value === '') return false
  if (positiveControl.value === negativeControl.value) return false
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
    positive_control: positiveControl.value,
    negative_control: negativeControl.value,
  })
}

const close = () => {
  emit('update:open', false)
}

const outputTitles = computed(() => ({
  completed: t('experiments.analysis.completed'),
  completedWithWarnings: t('experiments.analysis.completed_with_warnings'),
  failed: t('experiments.analysis.failed'),
}))

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

        <div class="grid grid-cols-2 gap-3">
          <div>
            <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
              {{ t('experiments.analysis.positive_control') }}
            </label>
            <select
              v-model="positiveControl"
              :class="[fieldClass, 'cursor-pointer']"
              :disabled="analysisStore.isRunning"
            >
              <option v-for="wellType in wellTypes" :key="`positive-${wellType}`" :value="wellType">
                {{ wellType }}
              </option>
            </select>
          </div>
          <div>
            <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
              {{ t('experiments.analysis.negative_control') }}
            </label>
            <select
              v-model="negativeControl"
              :class="[fieldClass, 'cursor-pointer']"
              :disabled="analysisStore.isRunning"
            >
              <option v-for="wellType in wellTypes" :key="`negative-${wellType}`" :value="wellType">
                {{ wellType }}
              </option>
            </select>
          </div>
        </div>
        <p v-if="positiveControl !== '' && positiveControl === negativeControl" class="pl-1 text-sm text-red-600">
          {{ t('experiments.analysis.same_controls') }}
        </p>

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

        <CommandOutputLog
          v-if="isThisExperiment && (analysisStore.messages.length > 0 || analysisStore.isRunning)"
          :messages="analysisStore.messages"
          :status="analysisStore.status"
          :titles="outputTitles"
        >
          <p v-if="analysisStore.isRunning" class="px-2 py-0.5 text-xs text-slate-500">
            {{ t('experiments.analysis.running') }}
          </p>
        </CommandOutputLog>
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
