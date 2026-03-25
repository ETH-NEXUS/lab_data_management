<script setup lang="ts">
import { getTextFromMessage } from '@nuxt/ui/utils/ai'
import type { ChatMessage, ChatMessageMetadata } from '~~/types/ml'
import { renderMarkdown } from '~~/utils/markdown'

const props = defineProps<{
  messages: ChatMessage[]
  status?: 'submitted' | 'streaming' | 'ready' | 'error'
}>()

const getDurationMs = (message: ChatMessage): number | null => {
  const metadata = message.metadata as ChatMessageMetadata | undefined
  return typeof metadata?.durationMs === 'number' ? metadata.durationMs : null
}

const getThinking = (message: ChatMessage): string | null => {
  const metadata = message.metadata as ChatMessageMetadata | undefined
  const thinking = metadata?.thinking
  if (typeof thinking !== 'string') return null
  return thinking.trim() || null
}

const formatDuration = (ms: number): string => {
  const s = ms / 1000
  return s < 10 ? `${s.toFixed(2)}s` : `${s.toFixed(1)}s`
}

const isEmptyAssistant = (message: ChatMessage): boolean => {
  if (message.role !== 'assistant') return false
  const parts = message.parts ?? []
  const text = parts
    .map((part) => part.text)
    .join('')
    .trim()
  return text.length === 0
}

const getRenderedContent = (message: ChatMessage): string => {
  const text = getTextFromMessage(message)
  return renderMarkdown(text)
}

const hasThinking = (message: ChatMessage): boolean => getThinking(message) != null
const hasAnyThinking = computed(() => props.messages.some((message) => hasThinking(message)))

const assistantProps = computed(() =>
  hasAnyThinking.value
    ? ({ variant: 'soft', icon: 'i-lucide-bot', ui: { content: 'w-full' } } as const)
    : ({ variant: 'soft', icon: 'i-lucide-bot' } as const),
)

const openThinkingIds = ref(new Set<string>())

const syncThinkingState = () => {
  const next = new Set(openThinkingIds.value)
  const lastAssistant = [...props.messages].reverse().find((message) => message.role === 'assistant')

  if (!lastAssistant) {
    openThinkingIds.value = next
    return
  }

  const lastMessageHasThinking = getThinking(lastAssistant) != null

  if (!lastMessageHasThinking) {
    openThinkingIds.value = next
    return
  }

  if (isEmptyAssistant(lastAssistant)) {
    next.add(lastAssistant.id)
  } else {
    next.delete(lastAssistant.id)
  }

  openThinkingIds.value = next
}

watch(
  () => props.messages,
  () => syncThinkingState(),
  { deep: true, immediate: true },
)

const isThinkingOpen = (messageId: string) => openThinkingIds.value.has(messageId)

const handleThinkingToggle = (messageId: string, values: string[]) => {
  const next = new Set(openThinkingIds.value)
  if (values.includes('thinking')) {
    next.add(messageId)
  } else {
    next.delete(messageId)
  }
  openThinkingIds.value = next
}
</script>

<template>
  <UChatMessages
    :messages="props.messages"
    :spacing-offset="80"
    :status="props.status"
    :user="{ variant: 'soft', icon: 'i-lucide-user' }"
    :assistant="assistantProps"
    compact
    class="h-full w-full"
  >
    <template #content="{ message }">
      <div class="flex flex-col gap-1">
        <UAccordion
          v-if="getThinking(message)"
          :items="
            [{ label: 'Thinking', value: 'thinking', content: getThinking(message) ?? '' }] as Array<{
              label: string
              value: string
              content: string
            }>
          "
          :model-value="isThinkingOpen(message.id) ? 'thinking' : undefined"
          size="xs"
          variant="ghost"
          class="w-max max-w-full text-[11px]"
          :unmount-on-hide="false"
          @update:model-value="(value) => handleThinkingToggle(message.id, value ? [value as string] : [])"
        >
          <template #content="{ item }">
            <div class="text-muted whitespace-pre-wrap">{{ (item as { content: string }).content }}</div>
          </template>
        </UAccordion>
        <div
          v-if="message.role === 'assistant' && !isEmptyAssistant(message)"
          class="markdown-content prose prose-sm dark:prose-invert max-w-none"
          v-html="getRenderedContent(message)"
        ></div>
        <div v-else-if="message.role === 'user'" class="whitespace-pre-wrap">
          {{ getTextFromMessage(message) }}
        </div>
        <div v-if="getDurationMs(message) != null" class="text-muted text-[11px]">
          {{ formatDuration(getDurationMs(message)!) }}
        </div>
      </div>
    </template>
  </UChatMessages>
</template>
