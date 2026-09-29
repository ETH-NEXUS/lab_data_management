<script setup lang="ts">
import { computed } from 'vue'
import type { CommandMessage, CommandStatus } from '~/types/management'
import { summarizeCommandErrors } from '~/utils/commandErrors'

/**
 * The output of a command that runs in the celery container (an import of the
 * management page, or a statistical analysis): a summary once it has ended,
 * with its errors, and all its messages in colors.
 *
 * Props example:
 * - `messages = [{ level: 'info', text: 'Processing...' }, { level: 'error', text: 'File not found' }]`
 * - `status = 'failed'`
 * - `titles = { completed: 'Command completed.', completedWithWarnings: '...', failed: '...' }`
 *
 * The default slot is shown below the messages, e.g. a "still running" line.
 */
const props = defineProps<{
  messages: CommandMessage[]
  status: CommandStatus | null
  titles: {
    completed: string
    completedWithWarnings: string
    failed: string
  }
}>()

const { t } = useI18n()

/**
 * The summary above the output, once the command has ended.
 *
 * Returned data example:
 * - `{ color: 'error', icon: 'i-lucide-circle-x', title: 'Command failed. The errors are marked red below.' }`
 */
const summary = computed(() => {
  const hasWarnings = props.messages.some((message) => message.level === 'warning')

  if (props.status === 'failed') {
    return { color: 'error' as const, icon: 'i-lucide-circle-x', title: props.titles.failed }
  }
  if (props.status === 'completed' && hasWarnings) {
    return { color: 'warning' as const, icon: 'i-lucide-triangle-alert', title: props.titles.completedWithWarnings }
  }
  if (props.status === 'completed') {
    return { color: 'success' as const, icon: 'i-lucide-circle-check', title: props.titles.completed }
  }
  return null
})

// The error texts, shown in the summary so they are seen without scrolling the log
const errors = computed(() => summarizeCommandErrors(props.messages))

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
  <div class="space-y-2">
    <UAlert v-if="summary" :color="summary.color" :icon="summary.icon" :title="summary.title" variant="subtle">
      <template v-if="props.status === 'failed' && errors.shown.length > 0" #description>
        <ul class="list-disc space-y-1 pl-5">
          <li v-for="(text, index) in errors.shown" :key="index" class="whitespace-pre-wrap" v-text="text" />
        </ul>
        <p v-if="errors.notShown > 0" class="mt-1">
          {{ t('common.command_output.more_errors', { count: errors.notShown }) }}
        </p>
      </template>
    </UAlert>
    <div class="max-h-80 overflow-auto rounded-xl border border-slate-200 bg-white p-2">
      <p
        v-for="(message, index) in props.messages"
        :key="index"
        class="rounded px-2 py-0.5 font-mono text-xs whitespace-pre-wrap"
        :class="messageClass(message)"
        v-text="message.text"
      />
      <slot />
    </div>
  </div>
</template>
