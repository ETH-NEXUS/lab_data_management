<script setup lang="ts">
import { storeToRefs } from 'pinia'
import { useMlChatStore } from '~~/stores/ml/chat'
import { useMlModelsStore } from '~~/stores/ml/models'

const chatStore = useMlChatStore()
const modelsStore = useMlModelsStore()

const { input, status, messages, durationLabel, requestStartedAtMs } = storeToRefs(chatStore)
const { models, selectedModel, selectedModelItem } = storeToRefs(modelsStore)

const scrollContainer = ref<HTMLElement | null>(null)

onMounted(() => {
  modelsStore.loadModels()
})

onBeforeUnmount(() => {
  chatStore.abortActiveRequest()
})

const scrollToBottom = async () => {
  await nextTick()
  const container = scrollContainer.value
  if (!container) return
  container.scrollTop = container.scrollHeight
}

watch(
  () => [messages.value, status.value],
  () => {
    scrollToBottom()
  },
  { deep: true },
)
</script>

<template>
  <div class="relative flex h-full w-full flex-col gap-2 sm:gap-3 sm:p-3">
    <div ref="scrollContainer" class="h-full min-h-0 flex-1 overflow-y-scroll pb-64">
      <MessageList :messages="messages" :status="status" />
    </div>

    <UChatPrompt
      v-model="input"
      variant="soft"
      placeholder="Write your question here ..."
      class="right-0 bottom-0 left-0 w-full text-sm sm:text-base"
      :ui="{ base: 'max-h-16 sm:max-h-20' }"
      @submit="chatStore.submit"
    >
      <UChatPromptSubmit :status="status" class="rounded-full" />

      <template #footer>
        <div class="flex w-full items-center justify-between gap-3">
          <USelect
            v-model="selectedModel"
            :items="models"
            :icon="selectedModelItem.icon"
            placeholder="Select a model"
            variant="ghost"
            size="sm"
            class="min-w-0 flex-1"
          />
          <span v-if="durationLabel" class="text-muted shrink-0 text-xs">
            {{ requestStartedAtMs != null ? 'Query' : 'Last' }}: {{ durationLabel }}
          </span>
        </div>
      </template>
    </UChatPrompt>
  </div>
</template>
