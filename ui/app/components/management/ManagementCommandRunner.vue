<script setup lang="ts">
import { onMounted } from 'vue'
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

onMounted(() => {
  managementStore.clearCommandOutput()
})

/**
 * Builds command payload with legacy-compatible management fields.
 *
 * Accepted data example:
 * - `formData = { input_file: '/data/file.csv', project_name: 'P1' }`
 *
 * Returned data example:
 * - `{ input_file: '/data/file.csv', project_name: 'P1', room_name: '12', command: 'import', what: 'library_plate', is_control_plate: true }`
 */
const buildCommandPayload = (formData: GeneralFormData): GeneralFormData => {
  const payload: GeneralFormData = { ...formData }

  const currentUserId = authStore.user?.id
  payload.room_name = currentUserId ? String(currentUserId) : `room_${Date.now()}`
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
  managementStore.commandOutput = `Executing command: ${props.command}\n...`
  const payload = buildCommandPayload(formData)
  await managementStore.runCommand(payload)
}
</script>

<template>
  <div class="space-y-4">
    <ManagementDynamicForm
      :options="props.options"
      :is-submitting="managementStore.isRunningCommand"
      @submit="onSubmit"
    />

    <div v-if="managementStore.commandOutput !== ''" class="space-y-2">
      <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
        {{ t('management.logs') }}
      </p>
      <div class="max-h-60 overflow-auto rounded-xl bg-slate-900 p-3 text-slate-100 shadow-inner">
        <pre class="font-mono text-xs whitespace-pre-wrap">{{ managementStore.commandOutput }}</pre>
      </div>
    </div>
  </div>
</template>
