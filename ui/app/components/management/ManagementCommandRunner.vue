<script setup lang="ts">
import { computed } from 'vue'
import CommandOutputLog from '~/components/common/CommandOutputLog.vue'
import ManagementDynamicForm from '~/components/management/ManagementDynamicForm.vue'
import { useAuthStore } from '~/stores/auth'
import { useManagementStore } from '~/stores/management'
import type { GeneralFormData, Options } from '~/types/lab'

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

const outputTitles = computed(() => ({
  completed: t('management.command_completed'),
  completedWithWarnings: t('management.command_completed_with_warnings'),
  failed: t('management.command_failed'),
}))
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
      <CommandOutputLog
        :messages="managementStore.commandMessages"
        :status="managementStore.commandStatus"
        :titles="outputTitles"
      />
    </div>
  </div>
</template>
