<script setup lang="ts">
import { computed } from 'vue'
import type { Experiment } from '~/types/experiments'
import type { Project } from '~/types/projects'

type UsageUnitOption = {
  id: number
  label: string
}

type Props = {
  open: boolean
  isSavingUsage: boolean
  stockItemLabel: string
  selectedProjectId: string
  selectedExperimentId: string
  quantityUsed: string
  selectedUsageUnitId: string
  usageNotes: string
  projects: Project[]
  filteredExperiments: Experiment[]
  usageUnitOptions: UsageUnitOption[]
}

const props = defineProps<Props>()

const emit = defineEmits<{
  (
    e:
      | 'update:selected-project-id'
      | 'update:selected-experiment-id'
      | 'update:quantity-used'
      | 'update:selected-usage-unit-id'
      | 'update:usage-notes',
    value: string,
  ): void
  (e: 'cancel' | 'save'): void
}>()

const { t } = useI18n()

const selectedProjectIdModel = computed<string>({
  get: () => props.selectedProjectId,
  set: (value) => emit('update:selected-project-id', value),
})

const selectedExperimentIdModel = computed<string>({
  get: () => props.selectedExperimentId,
  set: (value) => emit('update:selected-experiment-id', value),
})

const quantityUsedModel = computed<string>({
  get: () => props.quantityUsed,
  set: (value) => emit('update:quantity-used', value),
})

const selectedUsageUnitIdModel = computed<string>({
  get: () => props.selectedUsageUnitId,
  set: (value) => emit('update:selected-usage-unit-id', value),
})

const usageNotesModel = computed<string>({
  get: () => props.usageNotes,
  set: (value) => emit('update:usage-notes', value),
})
</script>

<template>
  <div v-if="props.open" class="space-y-3 rounded-lg border border-slate-200 bg-slate-50 p-3">
    <p class="text-xs text-slate-600">
      {{ t('inventory.stock_drawer.record_usage.stock_item', { label: props.stockItemLabel }) }}
    </p>

    <div class="grid gap-3 sm:grid-cols-2">
      <div class="space-y-1 sm:col-span-2">
        <label class="block text-sm font-medium text-slate-700">{{
          t('inventory.stock_drawer.record_usage.project')
        }}</label>
        <select
          v-model="selectedProjectIdModel"
          class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
          :disabled="props.isSavingUsage"
        >
          <option value="">{{ t('inventory.stock_drawer.record_usage.select_project') }}</option>
          <option v-for="project in props.projects" :key="project.id" :value="String(project.id)">
            {{ project.name }}
          </option>
        </select>
      </div>

      <div class="space-y-1 sm:col-span-2">
        <label class="block text-sm font-medium text-slate-700">{{
          t('inventory.stock_drawer.record_usage.experiment_optional')
        }}</label>
        <select
          v-model="selectedExperimentIdModel"
          class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
          :disabled="props.isSavingUsage"
        >
          <option value="">{{ t('inventory.stock_drawer.record_usage.no_experiment') }}</option>
          <option v-for="experiment in props.filteredExperiments" :key="experiment.id" :value="String(experiment.id)">
            {{ experiment.name }}
          </option>
        </select>
      </div>

      <div class="space-y-1">
        <label class="block text-sm font-medium text-slate-700">{{
          t('inventory.stock_drawer.record_usage.quantity_used')
        }}</label>
        <input
          v-model="quantityUsedModel"
          type="text"
          class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
          :disabled="props.isSavingUsage"
        />
      </div>

      <div class="space-y-1">
        <label class="block text-sm font-medium text-slate-700">{{
          t('inventory.stock_drawer.record_usage.usage_unit')
        }}</label>
        <select
          v-model="selectedUsageUnitIdModel"
          class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
          :disabled="props.isSavingUsage"
        >
          <option value="">{{ t('inventory.stock_drawer.record_usage.select_usage_unit') }}</option>
          <option v-for="unitOption in props.usageUnitOptions" :key="unitOption.id" :value="String(unitOption.id)">
            {{ unitOption.label }}
          </option>
        </select>
      </div>

      <div class="space-y-1 sm:col-span-2">
        <label class="block text-sm font-medium text-slate-700">Notes</label>
        <textarea
          v-model="usageNotesModel"
          rows="3"
          class="w-full rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 outline-none focus:border-blue-300 focus:ring-2 focus:ring-blue-300/50"
          :disabled="props.isSavingUsage"
        />
      </div>
    </div>

    <div class="flex justify-end gap-2">
      <UButton
        variant="ghost"
        color="neutral"
        :label="t('common.actions.cancel')"
        :disabled="props.isSavingUsage"
        @click="emit('cancel')"
      />
      <UButton
        color="primary"
        :label="t('common.actions.save')"
        :loading="props.isSavingUsage"
        :disabled="props.isSavingUsage"
        @click="emit('save')"
      />
    </div>
  </div>
</template>
