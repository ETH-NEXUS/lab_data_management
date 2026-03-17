<script setup lang="ts">
import { onMounted, ref } from 'vue'
import { useAPI } from '~/composables/useAPI'
import MessagesIntroCard from '~/components/messages/MessagesIntroCard.vue'
import MessagesProblematicPlatesCard from '~/components/messages/MessagesProblematicPlatesCard.vue'
import MessagesThresholdCard from '~/components/messages/MessagesThresholdCard.vue'
import MessagesThresholdModal from '~/components/messages/MessagesThresholdModal.vue'
import type { RedFlagInfo, Threshold, ThresholdListResponse, ThresholdUpdatePayload } from '~/types/messages'
import {
  DEFAULT_THRESHOLD,
  RECALCULATE_STATUS_ENDPOINT,
  RED_FLAG_ENDPOINT,
  THRESHOLDS_ENDPOINT,
  normalizeThreshold,
} from '~/utils/messages'

const { t } = useI18n()
const toast = useToast()

const isLoading = ref(false)
const isThresholdModalOpen = ref(false)
const isUpdatingThreshold = ref(false)
const isRecalculatingStatus = ref(false)

const threshold = ref<Threshold>(DEFAULT_THRESHOLD)
const redFlagInfo = ref<RedFlagInfo>({})

/**
 * Loads grouped red-flag warnings.
 *
 * Returned data example:
 * - `{ "Library A": { "PLATE-001": ["A01", "B02"] } }`
 */
const getInfo = async () => {
  const { data, error } = await useAPI<RedFlagInfo>(RED_FLAG_ENDPOINT, { method: 'GET' })
  if (error.value || !data.value) {
    console.error(error.value)
    return
  }

  redFlagInfo.value = data.value
}

/**
 * Loads current threshold values and stores the first result row.
 *
 * Returned data example:
 * - `{ results: [{ id: 1, dmso: 80, amount: 2.5 }] }`
 */
const getThresholds = async () => {
  const { data, error } = await useAPI<ThresholdListResponse>(THRESHOLDS_ENDPOINT, { method: 'GET' })
  if (error.value || !data.value) {
    console.error(error.value)
    return
  }

  const firstThreshold = data.value.results?.[0] ?? DEFAULT_THRESHOLD
  threshold.value = normalizeThreshold(firstThreshold)
}

const initializePage = async () => {
  isLoading.value = true
  try {
    await Promise.all([getInfo(), getThresholds()])
  } finally {
    isLoading.value = false
  }
}

/**
 * Updates threshold values and refreshes threshold display.
 *
 * Accepted payload example:
 * - `{ dmso: 75, amount: 2.0 }`
 */
const updateThreshold = async (payload: ThresholdUpdatePayload) => {
  isUpdatingThreshold.value = true
  try {
    const { error } = await useAPI<unknown>(`${THRESHOLDS_ENDPOINT}${threshold.value.id}/`, {
      method: 'PATCH',
      body: payload,
    })

    if (error.value) {
      console.error(error.value)
      return
    }

    isThresholdModalOpen.value = false
    await getThresholds()
  } finally {
    isUpdatingThreshold.value = false
  }
}

const recalculateStatus = async () => {
  if (isRecalculatingStatus.value) return

  isRecalculatingStatus.value = true
  try {
    toast.add({
      title: t('messages_page.recalculate.in_progress'),
      description: t('messages_page.recalculate.in_progress_description'),
      color: 'info',
      duration: 3000,
    })

    const { error } = await useAPI<unknown>(RECALCULATE_STATUS_ENDPOINT, { method: 'GET' })
    if (error.value) {
      console.error(error.value)
      return
    }

    await getInfo()
  } finally {
    isRecalculatingStatus.value = false
  }
}

onMounted(async () => {
  await initializePage()
})
</script>

<template>
  <section class="space-y-5 p-6">
    <h1 class="text-center text-4xl font-semibold text-blue-700">
      {{ t('messages_page.title') }}
    </h1>

    <MessagesIntroCard />

    <MessagesProblematicPlatesCard :red-flag-info="redFlagInfo" :is-loading="isLoading" />

    <MessagesThresholdCard
      :threshold="threshold"
      :is-recalculating="isRecalculatingStatus"
      @edit-threshold="isThresholdModalOpen = true"
      @recalculate-status="recalculateStatus"
    />

    <MessagesThresholdModal
      v-model:open="isThresholdModalOpen"
      :threshold="threshold"
      :is-saving="isUpdatingThreshold"
      @submit="updateThreshold"
    />
  </section>
</template>
