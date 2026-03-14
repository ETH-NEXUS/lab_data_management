<script setup lang="ts">
import { computed } from 'vue'
import type { RedFlagInfo } from '~/types/messages'

const props = defineProps<{
  redFlagInfo: RedFlagInfo
  isLoading?: boolean
}>()

const { t } = useI18n()

/**
 * Sorted entries for top-level library groups.
 *
 * Returned data example:
 * - `[["Library A", { "PLATE-001": ["A01"] }], ["Library B", { "PLATE-002": ["B03"] }]]`
 */
const libraryEntries = computed(() => {
  return Object.entries(props.redFlagInfo).sort(([libraryA], [libraryB]) => libraryA.localeCompare(libraryB))
})

/**
 * Returns sorted plate entries for one library section.
 *
 * Returned data example:
 * - `[["PLATE-001", ["A01", "B04"]], ["PLATE-002", ["H12"]]]`
 */
const getPlateEntries = (platesByBarcode: Record<string, string[]>) => {
  return Object.entries(platesByBarcode).sort(([plateA], [plateB]) => plateA.localeCompare(plateB))
}
</script>

<template>
  <UCard
    class="mx-auto w-[80%]"
    :ui="{
      root: 'border border-white/40 bg-white/30 backdrop-blur-md shadow-sm divide-y divide-white/20',
    }"
  >
    <template #header>
      <p class="text-primary font-semibold">{{ t('messages_page.sections.problematic_plates.title') }}</p>
    </template>

    <p v-if="props.isLoading" class="text-sm text-slate-600">
      {{ t('messages_page.loading') }}
    </p>

    <p v-else-if="libraryEntries.length === 0" class="text-sm text-slate-600">
      {{ t('messages_page.sections.problematic_plates.empty') }}
    </p>

    <div v-else class="space-y-2">
      <details
        v-for="[libraryName, platesByBarcode] in libraryEntries"
        :key="libraryName"
        class="group rounded-lg border border-slate-200 bg-slate-100/80 p-3"
      >
        <summary class="flex cursor-pointer list-none items-center justify-between gap-3">
          <span class="text-primary truncate text-sm font-semibold">{{ libraryName }}</span>
          <UIcon
            name="i-heroicons-chevron-right"
            class="size-5 shrink-0 text-slate-500 transition-transform duration-200 group-open:rotate-90"
          />
        </summary>

        <div class="mt-2 space-y-2">
          <details
            v-for="[plateBarcode, wells] in getPlateEntries(platesByBarcode)"
            :key="`${libraryName}-${plateBarcode}`"
            class="group rounded-md border border-slate-200 bg-white p-2"
          >
            <summary class="flex cursor-pointer list-none items-center justify-between gap-3">
              <span class="text-secondary truncate text-sm">{{ plateBarcode }}</span>
              <span class="text-xs text-slate-500">{{ wells.length }}</span>
            </summary>

            <ul class="mt-2 space-y-1">
              <li
                v-for="well in wells"
                :key="`${libraryName}-${plateBarcode}-${well}`"
                class="rounded-md bg-slate-50 px-2 py-1 text-sm text-slate-700"
              >
                {{ well }}
              </li>
            </ul>
          </details>
        </div>
      </details>
    </div>
  </UCard>
</template>
