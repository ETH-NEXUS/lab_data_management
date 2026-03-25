<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import BaseField from '~/components/common/BaseField.vue'
import type { GeneralFormData, Options } from '~/types/lab'

type Props = {
  options: Options
  isSubmitting?: boolean
}

const props = withDefaults(defineProps<Props>(), {
  isSubmitting: false,
})

const emit = defineEmits<{
  (e: 'submit', formData: GeneralFormData): void
}>()

const toast = useToast()
const { t } = useI18n()
const formState = ref<GeneralFormData>({})

const optionEntries = computed(() => {
  return Object.entries(props.options)
})

/**
 * Resets reactive form values from one options map.
 *
 * Accepted data example:
 * - `{ machine: { type: 'str', required: true }, dry_run: { type: 'bool', required: false } }`
 *
 * Result example:
 * - `formState = { machine: '', dry_run: false }`
 */
const initializeForm = (): void => {
  const nextState: GeneralFormData = {}

  for (const [key, option] of optionEntries.value) {
    if (option.type === 'bool') {
      nextState[key] = false
    } else {
      nextState[key] = ''
    }
  }

  formState.value = nextState
}

watch(
  () => props.options,
  () => {
    initializeForm()
  },
  { immediate: true, deep: true },
)

/**
 * Validates required fields before submit.
 *
 * Returned data example:
 * - `{ isValid: false, firstInvalidLabel: 'Experiment name' }`
 * - `{ isValid: true, firstInvalidLabel: null }`
 */
const validateForm = (): { isValid: boolean; firstInvalidLabel: string | null } => {
  for (const [key, option] of optionEntries.value) {
    if (!option.required) {
      continue
    }

    const value = formState.value[key]

    if (option.type === 'bool') {
      if (!value) {
        return { isValid: false, firstInvalidLabel: option.label }
      }

      continue
    }

    if (typeof value !== 'string' || value.trim() === '') {
      return { isValid: false, firstInvalidLabel: option.label }
    }
  }

  return { isValid: true, firstInvalidLabel: null }
}

const submitForm = (): void => {
  const validation = validateForm()
  if (!validation.isValid) {
    toast.add({
      title: t('management.validation.required_title'),
      description: t('management.validation.required_field', { label: validation.firstInvalidLabel ?? '' }),
      color: 'error',
      duration: 3000,
    })
    return
  }

  emit('submit', { ...formState.value })
}

/**
 * Enables dropping plain text (filesystem path) onto string fields.
 *
 * Accepted data example:
 * - drag payload: `'/data/imports/file.csv'`
 * - target field key: `'input_file'`
 */
const onDragOverStringField = (event: DragEvent): void => {
  event.preventDefault()

  if (event.dataTransfer) {
    event.dataTransfer.dropEffect = 'copy'
  }
}

/**
 * Reads plain-text payload from one drag event.
 *
 * Returned value examples:
 * - `'/data/imports/file.csv'`
 * - `''` when drag payload has no text path
 */
const readDroppedPath = (event: DragEvent): string => {
  if (!event.dataTransfer) {
    return ''
  }

  const payloadCandidates = [
    event.dataTransfer.getData('application/x-management-path'),
    event.dataTransfer.getData('text/plain'),
    event.dataTransfer.getData('text'),
  ]

  for (const payload of payloadCandidates) {
    const normalizedPayload = payload.trim()
    if (normalizedPayload !== '') {
      return normalizedPayload
    }
  }

  return ''
}

/**
 * Writes dropped text payload into one form field.
 *
 * Accepted data examples:
 * - `key = 'path'`, payload = `'/data/mappings'`
 * - `key = 'input_file'`, payload = `'/data/imports/control.csv'`
 */
const onDropToStringField = (key: string, event: DragEvent): void => {
  event.preventDefault()

  const droppedValue = readDroppedPath(event)
  if (droppedValue === '') {
    return
  }

  formState.value[key] = droppedValue
}
</script>

<template>
  <div class="space-y-4">
    <div v-for="[key, option] in optionEntries" :key="`management-option-${key}`" class="space-y-2">
      <div v-if="option.type === 'str' && option.choices" class="space-y-2">
        <label class="mb-2 block pl-1 text-sm font-semibold tracking-[0.04em] text-slate-600">
          {{ option.label }}
        </label>

        <select
          v-model="formState[key]"
          class="w-full rounded-md border border-slate-200 bg-white px-6 py-4 text-slate-800 shadow-[0_1px_0_rgba(15,23,42,0.02)] ring-offset-0 transition outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
        >
          <option value="">{{ t('management.form.select_placeholder') }}</option>
          <option v-for="choice in option.choices" :key="`${key}-${choice}`" :value="choice">
            {{ choice }}
          </option>
        </select>
      </div>

      <div
        v-else-if="option.type === 'str'"
        class="rounded-lg"
        @dragenter="onDragOverStringField"
        @dragover="onDragOverStringField"
        @drop="onDropToStringField(key, $event)"
      >
        <BaseField
          v-model="formState[key]"
          :label="option.label"
          :placeholder="option.label"
          @dragenter="onDragOverStringField"
          @dragover="onDragOverStringField"
          @drop="onDropToStringField(key, $event)"
        />
      </div>

      <label v-else class="inline-flex cursor-pointer items-center gap-2">
        <input
          v-model="formState[key]"
          type="checkbox"
          class="h-4 w-4 rounded border-slate-300 text-blue-600 focus:ring-blue-300"
        />
        <span class="text-sm text-slate-700">{{ option.label }}</span>
      </label>
    </div>

    <BaseButton
      :label="t('action.submit')"
      :on-click="submitForm"
      variant="primary"
      size="sm"
      width="auto"
      :loading="props.isSubmitting"
      :disabled="props.isSubmitting"
    />
  </div>
</template>
