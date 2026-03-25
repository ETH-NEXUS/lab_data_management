<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import BaseField from '~/components/common/BaseField.vue'
import WavesModalWrapper from '~/components/common/WavesModalWrapper.vue'
import { useExperimentStore } from '~/stores/experiments'
import { getErrorMessage } from '~/utils/errors'

const props = defineProps<{
  open: boolean
  experimentId: number
  labels: string[]
  measurementTimestamps?: {
    [key: string]: string[]
  }
}>()

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
  (e: 'calculated'): void
}>()

const { t } = useI18n()
const toast = useToast()
const experimentStore = useExperimentStore()

const buttons: string[] = ['1', '2', '3', '+', '4', '5', '6', '-', '7', '8', '9', '*', 'C', '0', '/', 'log']
const currentExpression = ref('')
const previousButton = ref('')
const newLabel = ref('')
const usedLabels = ref<string[]>([])
const combinedLabels = ref<string[]>([])

watch(
  () => props.open,
  (isOpen) => {
    if (!isOpen) return

    currentExpression.value = ''
    previousButton.value = ''
    newLabel.value = ''
    usedLabels.value = []
    combinedLabels.value = createCombinedLabels()
  },
  { immediate: true },
)

/**
 * Checks whether one measurement label has multiple timestamps.
 *
 * Accepted input examples:
 * - `['2026-01-01T10:00:00Z', '2026-01-01T10:00:00Z']` -> `false`
 * - `['2026-01-01T10:00:00Z', '2026-01-01T11:00:00Z']` -> `true`
 */
const hasDifferentTimestamps = (timestamps: string[]) => {
  const set = new Set(timestamps)
  return set.size !== 1
}

/**
 * Builds labels in the legacy "Label --> Timestamp" format.
 *
 * Returned data example:
 * - `['OD600 --> 2026-01-01T10:00:00Z', 'OD600 --> 2026-01-01T11:00:00Z']`
 */
const createCombinedLabels = () => {
  const mergedLabels: string[] = []

  if (!props.measurementTimestamps) {
    return mergedLabels
  }

  for (const [label, timestamps] of Object.entries(props.measurementTimestamps)) {
    if (!hasDifferentTimestamps(timestamps)) {
      continue
    }

    for (const timestamp of timestamps) {
      // Keep this separator
      mergedLabels.push(`${label} --> ${timestamp}`)
    }
  }

  return mergedLabels
}

const handleButtonClick = (button: string) => {
  switch (button) {
    case 'C':
      currentExpression.value = ''
      break
    case 'log':
      currentExpression.value += 'log('
      previousButton.value = button
      break
    default:
      if (previousButton.value === 'log') {
        currentExpression.value += `${button})`
      } else {
        currentExpression.value += button
      }
      previousButton.value = button
      break
  }
}

const handleDragStart = (label: string, event: DragEvent) => {
  event.dataTransfer?.setData('text/plain', label)
}

const handleDrop = (event: DragEvent) => {
  const label = event.dataTransfer?.getData('text/plain')
  if (!label) return

  currentExpression.value += label
  if (previousButton.value === 'log') {
    currentExpression.value += ')'
  }
  usedLabels.value.push(label)
}

const close = () => emit('update:open', false)

/**
 * Validates old calculator rules and sends the calculate request.
 *
 * Sent request example:
 * - `addNewMeasurement(null, 'A+B', 'Combined', ['A', 'B'], 21)`
 */
const addMeasurement = async () => {
  if (!newLabel.value) {
    toast.add({
      title: t('experiments.measurement_calculator.error_enter_name'),
      color: 'error',
      duration: 3500,
    })
    return
  }

  if (props.labels.includes(newLabel.value)) {
    toast.add({
      title: t('experiments.measurement_calculator.error_use_different_name'),
      color: 'error',
      duration: 3500,
    })
    return
  }

  const usedLabelsSet = new Set(usedLabels.value)
  const labelsSet = new Set<string>()
  for (const label of props.labels) {
    labelsSet.add(label)
  }
  for (const label of combinedLabels.value) {
    labelsSet.add(label)
  }

  const intersection = new Set<string>()
  for (const usedLabel of usedLabelsSet) {
    if (labelsSet.has(usedLabel)) {
      intersection.add(usedLabel)
    }
  }

  if (intersection.size === 0) {
    toast.add({
      title: t('experiments.measurement_calculator.error_use_existing_measurement'),
      color: 'error',
      duration: 4000,
    })
    return
  }

  const hasTimeSeriesLabels = usedLabels.value.some((label) => label.includes('-->'))
  const hasNormalLabels = usedLabels.value.some((label) => !label.includes('-->'))

  if (hasTimeSeriesLabels && hasNormalLabels) {
    toast.add({
      title: t('experiments.measurement_calculator.error_mixed_series_selection'),
      color: 'error',
      duration: 5000,
    })
    return
  }

  try {
    await experimentStore.addNewMeasurement(
      null,
      currentExpression.value,
      newLabel.value,
      usedLabels.value,
      props.experimentId,
    )

    emit('calculated')
    close()

    toast.add({
      title: t('experiments.measurement_calculator.success_saved'),
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: t('experiments.measurement_calculator.error_calculation_failed'),
      description: getErrorMessage(err),
      color: 'error',
      duration: 4500,
    })
  }
}

const canApply = computed(() => {
  if (experimentStore.isAddingMeasurement) return false
  if (newLabel.value.trim() === '') return false
  if (currentExpression.value.trim() === '') return false
  return true
})
</script>

<template>
  <WavesModalWrapper
    :open="props.open"
    :title="t('experiments.measurement_calculator.title')"
    :description="t('experiments.measurement_calculator.description')"
    :dismissible="!experimentStore.isAddingMeasurement"
    modal-class="w-full sm:max-w-4xl"
    body-container-class="w-full max-w-3xl px-8 pt-10"
    @update:open="emit('update:open', $event)"
  >
    <template #body>
      <div>
        <div class="calculator">
          <input
            v-model="currentExpression"
            type="text"
            class="calculator__input"
            @drop.prevent="handleDrop($event)"
            @dragover.prevent
          />
          <div class="calculator__buttons">
            <button
              v-for="button in buttons"
              :key="button"
              type="button"
              class="calculator__button"
              @click="handleButtonClick(button)"
            >
              {{ button }}
            </button>
          </div>
        </div>

        <div v-if="combinedLabels.length > 0" class="centered text-body2 text-blue-600">
          {{ t('experiments.measurement_calculator.labels_all_series') }}
        </div>
        <div class="labels">
          <div
            v-for="(label, index) in props.labels"
            :key="`${label}-${index}`"
            class="label"
            draggable="true"
            @dragstart="handleDragStart(label, $event)"
          >
            {{ label }}
          </div>
        </div>

        <div v-if="combinedLabels.length > 0" class="centered text-body2 text-blue-600">
          {{ t('experiments.measurement_calculator.labels_separate_series') }}
        </div>
        <div v-if="combinedLabels.length > 0" class="labels">
          <div
            v-for="(label, index) in combinedLabels"
            :key="`${label}-${index}`"
            class="label"
            draggable="true"
            @dragstart="handleDragStart(label, $event)"
          >
            {{ label }}
          </div>
        </div>

        <div class="mt-6">
          <BaseField
            v-model="newLabel"
            :label="t('experiments.measurement_calculator.new_measurement_name_label')"
            :placeholder="t('experiments.measurement_calculator.new_measurement_name_placeholder')"
          />
        </div>
      </div>
    </template>

    <template #footer>
      <BaseButton
        :label="t('common.actions.cancel')"
        :on-click="close"
        variant="secondary"
        size="sm"
        width="auto"
        :disabled="experimentStore.isAddingMeasurement"
      />
      <BaseButton
        :label="t('experiments.measurement_calculator.apply_button')"
        :on-click="addMeasurement"
        variant="primary"
        size="sm"
        width="auto"
        :loading="experimentStore.isAddingMeasurement"
        :disabled="!canApply"
      />
    </template>
  </WavesModalWrapper>
</template>

<style scoped>
.centered {
  text-align: center;
  margin: 10px 0;
}

.calculator {
  display: flex;
  flex-direction: column;
  align-items: center;
  padding: 20px;
  margin: 20px 0;
  background-color: #f2f2f2;
  border-radius: 10px;
}

.calculator__input {
  margin-bottom: 10px;
  padding: 10px;
  font-size: 14px;
  text-align: right;
  width: 100%;
  box-sizing: border-box;
  background: #ffffff;
}

.calculator__buttons {
  display: grid;
  grid-template-columns: repeat(4, 1fr);
  gap: 10px;
}

.calculator__button {
  background-color: #ffffff;
  color: #000000;
  border: none;
  border-radius: 5px;
  padding: 10px;
  font-size: 18px;
  cursor: pointer;
  transition: all 0.2s ease-in-out;
}

.calculator__button:hover {
  background-color: #cccccc;
}

.labels {
  display: flex;
  flex-wrap: wrap;
  justify-content: center;
  align-items: center;
  margin-top: 20px;
}

.label {
  padding: 10px 20px;
  margin: 10px;
  border-radius: 20px;
  background-color: #f2f2f2;
  font-size: 12px;
  font-weight: bold;
  color: #333333;
  box-shadow: 2px 2px 4px rgb(0 0 0 / 20%);
  cursor: move;
}

.label:hover {
  background-color: #cccccc;
}

.label:active {
  box-shadow: none;
  background-color: rgb(255 255 255 / 50%);
}
</style>
