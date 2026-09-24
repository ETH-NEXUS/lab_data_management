<script setup lang="ts">
import { computed } from 'vue'
import ManagementDynamicForm from '~/components/management/ManagementDynamicForm.vue'
import { useAuthStore } from '~/stores/auth'
import { useManagementStore } from '~/stores/management'
import type { GeneralFormData, Options } from '~/types/lab'
import type { CommandMessage } from '~/types/management'
import { summarizeCommandErrors } from '~/utils/commandErrors'

type Props = {
  options: Options
  command: string
  what: string
}

const props = defineProps<Props>()

const managementStore = useManagementStore()
const authStore = useAuthStore()
const { t } = useI18n()

/**
 * Builds command payload with legacy-compatible management fields.
 *
 * Accepted data example:
 * - `formData = { input_file: '/data/file.csv', project_name: 'P1' }`
 *
 * Returned data example:
 * - `{ input_file: '/data/file.csv', project_name: 'P1', room_name: '12_1726563600000', command: 'import', what: 'library_plate', is_control_plate: true }`
 */
const buildCommandPayload = (formData: GeneralFormData): GeneralFormData => {
  const payload: GeneralFormData = { ...formData }

  const currentUserId = authStore.user?.id
  // A new room for every run, so the output of an earlier run is never shown
  payload.room_name = `${currentUserId ?? 'room'}_${Date.now()}`
  payload.command = props.command

  if (props.what && props.command === 'import') {
    if (props.what === 'control_plate') {
      payload.what = 'library_plate'
      payload.is_control_plate = true
    } else if (props.what === 'library_plate') {
      payload.what = 'library_plate'
      payload.is_control_plate = false
    } else {
      payload.what = props.what
    }
  }

  return payload
}

const onSubmit = async (formData: GeneralFormData): Promise<void> => {
  const payload = buildCommandPayload(formData)
  await managementStore.runCommand(payload)
}

/**
 * The summary above the output, once the command has ended.
 *
 * Returned data example:
 * - `{ color: 'error', icon: 'i-lucide-circle-x', title: 'Command failed. The errors are marked red below.' }`
 */
const commandSummary = computed(() => {
  const hasWarnings = managementStore.commandMessages.some((message) => message.level === 'warning')

  if (managementStore.commandStatus === 'failed') {
    return { color: 'error' as const, icon: 'i-lucide-circle-x', title: t('management.command_failed') }
  }
  if (managementStore.commandStatus === 'completed' && hasWarnings) {
    return {
      color: 'warning' as const,
      icon: 'i-lucide-triangle-alert',
      title: t('management.command_completed_with_warnings'),
    }
  }
  if (managementStore.commandStatus === 'completed') {
    return { color: 'success' as const, icon: 'i-lucide-circle-check', title: t('management.command_completed') }
  }
  return null
})

// The error texts, shown in the summary so they are seen without scrolling the log
const commandErrors = computed(() => summarizeCommandErrors(managementStore.commandMessages))

const messageClass = (message: CommandMessage): string => {
  if (message.level === 'error') {
    return 'bg-red-50 text-red-800'
  }
  if (message.level === 'warning') {
    return 'bg-amber-50 text-amber-800'
  }
  if (message.level === 'success') {
    return 'text-green-700'
  }
  return 'text-slate-700'
}
</script>

<template>
  <div class="space-y-4">
    <ManagementDynamicForm
      :options="props.options"
      :is-submitting="managementStore.isRunningCommand"
      @submit="onSubmit"
    />

    <div v-if="managementStore.commandMessages.length > 0" class="space-y-2">
      <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
        {{ t('management.logs') }}
      </p>
      <UAlert
        v-if="commandSummary"
        :color="commandSummary.color"
        :icon="commandSummary.icon"
        :title="commandSummary.title"
        variant="subtle"
      >
        <template v-if="managementStore.commandStatus === 'failed' && commandErrors.shown.length > 0" #description>
          <ul class="list-disc space-y-1 pl-5">
            <li v-for="(text, index) in commandErrors.shown" :key="index" class="whitespace-pre-wrap" v-text="text" />
          </ul>
          <p v-if="commandErrors.notShown > 0" class="mt-1">
            {{ t('management.command_more_errors', { count: commandErrors.notShown }) }}
          </p>
        </template>
      </UAlert>
      <div class="max-h-80 overflow-auto rounded-xl border border-slate-200 bg-white p-2">
        <p
          v-for="(message, index) in managementStore.commandMessages"
          :key="index"
          class="rounded px-2 py-0.5 font-mono text-xs whitespace-pre-wrap"
          :class="messageClass(message)"
          v-text="message.text"
        />
      </div>
    </div>
  </div>
</template>
