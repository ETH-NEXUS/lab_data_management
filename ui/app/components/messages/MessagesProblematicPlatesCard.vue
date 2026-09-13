<script setup lang="ts">
import { computed } from 'vue'
import type { ProblematicWell, ProblematicWellReason, RedFlagInfo } from '~/types/messages'

const props = defineProps<{
  redFlagInfo: RedFlagInfo
  isLoading?: boolean
}>()

const { t } = useI18n()

/**
 * Sorted entries for top-level library groups.
 *
 * Returned data example:
 * - `[["Library A", { "PLATE-001": [{ position: "A01", ... }] }], ["Library B", { "PLATE-002": [] }]]`
 */
const libraryEntries = computed(() => {
  return Object.entries(props.redFlagInfo).sort(([libraryA], [libraryB]) => libraryA.localeCompare(libraryB))
})

/**
 * Returns sorted plate entries for one library section.
 *
 * Returned data example:
 * - `[["PLATE-001", [{ position: "A01", ... }]], ["PLATE-002", []]]`
 */
const getPlateEntries = (platesByBarcode: Record<string, ProblematicWell[]>) => {
  return Object.entries(platesByBarcode).sort(([plateA], [plateB]) => plateA.localeCompare(plateB))
}

/**
 * Formats the reported volume, or says that nothing was reported.
 *
 * Returned output examples:
 * - `formatVolume(1.31)` -> `'1.31 µL'`
 * - `formatVolume(null)` -> `'not reported'`
 */
const formatVolume = (value: number | null) => {
  if (value === null) return t('messages_page.sections.problematic_plates.not_reported')
  return `${value} ${t('unit.mikro')}`
}

/**
 * Formats the reported DMSO share, or says that nothing was reported.
 *
 * Returned output examples:
 * - `formatDmso(94.5)` -> `'94.5%'`
 * - `formatDmso(null)` -> `'not reported'`
 */
const formatDmso = (value: number | null) => {
  if (value === null) return t('messages_page.sections.problematic_plates.not_reported')
  return `${value}%`
}

/**
 * Names the thresholds a well is below, in words.
 *
 * Returned output example:
 * - `describeReasons({ reasons: ['volume', 'dmso'], ... })` -> `'volume below threshold, DMSO below threshold'`
 */
const describeReasons = (well: ProblematicWell) => {
  const descriptions = well.reasons.map((reason) => t(`messages_page.sections.problematic_plates.reasons.${reason}`))
  return descriptions.join(', ')
}

/**
 * Tells whether a well was marked because of the given threshold.
 * Used to highlight the value that is below it.
 */
const isBelow = (well: ProblematicWell, reason: ProblematicWellReason) => {
  return well.reasons.includes(reason)
}
</script>

<template>
  <UCard
    class="mx-auto w-[80%]"
    :ui="{
      root: 'core-card divide-y divide-slate-200/70',
    }"
  >
    <template #header>
      <p class="font-semibold text-blue-700">{{ t('messages_page.sections.problematic_plates.title') }}</p>
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
          <span class="truncate text-sm font-semibold text-blue-700">{{ libraryName }}</span>
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
                :key="`${libraryName}-${plateBarcode}-${well.position}`"
                class="flex flex-wrap items-baseline gap-x-3 rounded-md bg-slate-50 px-2 py-1 text-sm text-slate-700"
              >
                <span class="font-mono font-semibold">{{ well.position }}</span>
                <!-- The value below the threshold is shown in red and also named in words,
                     so the reason does not depend on seeing the color. -->
                <span :class="isBelow(well, 'volume') ? 'font-semibold text-red-700' : ''">
                  {{ t('messages_page.sections.problematic_plates.volume') }}:
                  {{ formatVolume(well.current_amount) }}
                </span>
                <span :class="isBelow(well, 'dmso') ? 'font-semibold text-red-700' : ''">
                  {{ t('messages_page.sections.problematic_plates.dmso') }}:
                  {{ formatDmso(well.current_dmso) }}
                </span>
                <span v-if="well.reasons.length > 0" class="text-xs text-slate-500">
                  {{ describeReasons(well) }}
                </span>
              </li>
            </ul>
          </details>
        </div>
      </details>
    </div>
  </UCard>
</template>
