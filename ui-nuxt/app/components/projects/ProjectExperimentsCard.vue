<script setup lang="ts">
import { computed, ref } from 'vue'
import BaseButton from '~/components/common/BaseButton.vue'
import { useProjectStore } from '~/stores/projects'
import type { Experiment } from '~/types/experiments'
import type { Project } from '~/types/projects'
import { getErrorMessage } from '~/utils/errors'

const props = defineProps<{
  project: Project
}>()

const emit = defineEmits<{
  (e: 'moved'): void
}>()

const { t } = useI18n()
const toast = useToast()
const projectStore = useProjectStore()

const selectedPlates = ref<string[]>([])
const selectedExperiment = ref<string | null>(null)

/**
 * Local project experiments list.
 *
 * Returned data examples:
 * - `[{ id: 3, name: 'Dose response', plates: [{ id: 12, barcode: 'A001' }] }]`
 * - `[]`
 */
const experiments = computed<Experiment[]>(() => {
  const experimentList: Experiment[] = []

  if (!props.project.experiments) {
    return experimentList
  }

  for (const experiment of props.project.experiments) {
    experimentList.push(experiment)
  }

  return experimentList
})

/**
 * Select options for the target experiment dropdown.
 *
 * Returned data examples:
 * - `['Dose response', 'QC run']`
 * - `[]`
 */
const experimentOptions = computed<string[]>(() => {
  const options: string[] = []

  for (const experiment of experiments.value) {
    options.push(experiment.name)
  }

  return options
})

/**
 * Disable "Move plates" button when one required input is missing.
 *
 * Returned data examples:
 * - `true` (nothing selected)
 * - `false` (target + at least one plate selected)
 */
const isMoveDisabled = computed<boolean>(() => {
  if (projectStore.isMovingPlates) return true
  if (!selectedExperiment.value) return true
  if (selectedPlates.value.length === 0) return true
  return false
})

const resetForm = () => {
  selectedPlates.value = []
  selectedExperiment.value = null
}

/**
 * Sends selected plate barcodes to backend and assigns them to one experiment name.
 *
 * Sent payload examples:
 * - `{ plate_barcodes: ['A001', 'A002'], experiment: 'Dose response' }`
 * - `{ plate_barcodes: ['CTRL-01'], experiment: 'QC run' }`
 */
const movePlates = async () => {
  if (selectedPlates.value.length === 0) {
    toast.add({
      title: t('projects.experiments_card.select_plate_warning'),
      color: 'warning',
      duration: 2500,
    })
    return
  }

  if (!selectedExperiment.value) {
    toast.add({
      title: t('projects.experiments_card.select_experiment_warning'),
      color: 'warning',
      duration: 2500,
    })
    return
  }

  try {
    await projectStore.movePlatesToExperiment({
      plate_barcodes: selectedPlates.value,
      experiment: selectedExperiment.value,
    })

    resetForm()
    emit('moved')

    toast.add({
      title: t('projects.experiments_card.move_success'),
      color: 'success',
      duration: 2500,
    })
  } catch (err: unknown) {
    toast.add({
      title: t('projects.experiments_card.move_failed'),
      description: getErrorMessage(err),
      color: 'error',
      duration: 4000,
    })
  }
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
      <div class="space-y-1">
        <p class="font-semibold">{{ t('projects.experiments_card.title') }}</p>
        <p class="text-xs text-slate-600">{{ t('projects.experiments_card.description') }}</p>
      </div>
    </template>

    <p v-if="experiments.length === 0" class="text-sm text-slate-600">
      {{ t('projects.experiments_card.no_experiments') }}
    </p>

    <section v-else class="flex flex-wrap gap-6">
      <article class="w-full max-w-[520px] space-y-2">
        <details
          v-for="experiment in experiments"
          :key="`experiment-${experiment.id}`"
          class="group rounded-lg border border-slate-200 bg-white/90 p-3"
        >
          <summary class="flex cursor-pointer list-none items-center justify-between gap-3">
            <span class="truncate text-sm font-semibold text-slate-800">{{ experiment.name }}</span>
            <span class="rounded bg-slate-100 px-2 py-0.5 text-xs text-slate-600">
              {{ experiment.plates.length }}
            </span>
          </summary>

          <div class="mt-2 space-y-2">
            <p v-if="experiment.plates.length === 0" class="text-xs text-slate-500">
              {{ t('projects.experiments_card.no_plates') }}
            </p>

            <label
              v-for="plate in experiment.plates"
              :key="`plate-${plate.id}`"
              class="flex cursor-pointer items-center gap-2 rounded-md bg-slate-50 px-2 py-1"
            >
              <input v-model="selectedPlates" type="checkbox" :value="plate.barcode" class="h-4 w-4" />
              <span class="font-mono text-sm text-slate-700">{{ plate.barcode }}</span>
            </label>
          </div>
        </details>
      </article>

      <article class="w-full max-w-[360px] space-y-3">
        <label class="block text-sm font-medium text-slate-700">
          {{ t('projects.experiments_card.choose_experiment') }}
        </label>

        <select
          v-model="selectedExperiment"
          class="w-full rounded-md border border-slate-300 bg-white px-3 py-2 text-sm text-slate-700"
        >
          <option :value="null">{{ t('projects.experiments_card.select_experiment_placeholder') }}</option>
          <option v-for="option in experimentOptions" :key="`option-${option}`" :value="option">
            {{ option }}
          </option>
        </select>

        <p class="text-xs text-slate-600">
          {{ t('projects.experiments_card.selected_plates_count', { count: selectedPlates.length }) }}
        </p>

        <BaseButton
          :label="t('projects.experiments_card.move_button')"
          :on-click="movePlates"
          variant="outline"
          size="sm"
          width="auto"
          class-name="text-teal-800 hover:text-teal-900 hover:bg-teal-50 hover:border-teal-200"
          :loading="projectStore.isMovingPlates"
          :disabled="isMoveDisabled"
        />
      </article>
    </section>
  </UCard>
</template>
