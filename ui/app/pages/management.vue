<script setup lang="ts">
import { computed, onMounted, ref } from 'vue'
import ManagementCommandRunner from '~/components/management/ManagementCommandRunner.vue'
import ManagementSelectedPaths from '~/components/management/ManagementSelectedPaths.vue'
import { useManagementStore } from '~/stores/management'
import type { ManagementTabKey } from '~/utils/managementCommandOptions'
import { managementCommandConfigs } from '~/utils/managementCommandOptions'
import { getErrorMessage } from '~/utils/errors'

const { t } = useI18n()
const managementStore = useManagementStore()

const activeTab = ref<ManagementTabKey>('map')
const pageErrorMessage = ref<string | null>(null)

const activeConfig = computed(() => {
  for (const config of managementCommandConfigs) {
    if (config.tab === activeTab.value) {
      return config
    }
  }

  return managementCommandConfigs[0]
})

onMounted(async () => {
  pageErrorMessage.value = null

  try {
    await managementStore.initialize()
  } catch (err: unknown) {
    pageErrorMessage.value = getErrorMessage(err)
  }
})

/**
 * Builds translated tab label from one command config.
 *
 * Accepted data example:
 * - `{ labelKey: 'management.echo', tab: 'map' }`
 *
 * Returned data example:
 * - `'Echo/M1000/Micro'`
 */
const getTabLabel = (labelKey: string): string => {
  return t(labelKey)
}
</script>

<template>
  <section class="space-y-6 p-6">
    <h1 class="text-center text-4xl font-semibold text-blue-700">
      {{ t('management.management') }}
    </h1>

    <p v-if="pageErrorMessage" class="rounded-md border border-red-200 bg-red-50 p-3 text-sm text-red-700">
      {{ pageErrorMessage }}
    </p>

    <ManagementSelectedPaths />

    <UCard
      class="mx-auto w-full max-w-5xl"
      :ui="{
        root: 'core-card divide-y divide-slate-200/70',
      }"
    >
      <template #header>
        <div class="flex flex-wrap gap-2">
          <button
            v-for="config in managementCommandConfigs"
            :key="`management-tab-${config.tab}`"
            type="button"
            class="inline-flex items-center justify-center rounded-full border px-4 py-2 text-xs font-semibold tracking-[0.08em] uppercase transition"
            :class="
              config.tab === activeTab
                ? 'border-blue-700 bg-blue-700 text-white'
                : 'border-blue-100 bg-white text-blue-700 hover:border-blue-200 hover:bg-blue-50'
            "
            @click="activeTab = config.tab"
          >
            {{ getTabLabel(config.labelKey) }}
          </button>
        </div>
      </template>

      <ManagementCommandRunner
        v-if="activeConfig"
        :options="activeConfig.options"
        :command="activeConfig.command"
        :what="activeConfig.what"
      />
    </UCard>
  </section>
</template>
