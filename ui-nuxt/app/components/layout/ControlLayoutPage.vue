<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import { useQueryClient } from '@tanstack/vue-query'
import BaseButton from '~/components/common/BaseButton.vue'
import BaseField from '~/components/common/BaseField.vue'
import WavesModalWrapper from '~/components/common/WavesModalWrapper.vue'
import DynamicPlate from '~/components/plates/DynamicPlate.vue'
import { useProjectQuery } from '~/composables/useProjectsQuery'
import { useLayoutStore } from '~/stores/layout'
import { PROJECTS_QUERY_KEY } from '~/types/projects'

const route = useRoute()
const toast = useToast()
const queryClient = useQueryClient()
const layoutStore = useLayoutStore()

const addLayoutDialog = ref(false)
const newBarcode = ref('')
const oldBarcode = ref('')

const projectId = computed<number | null>(() => {
  const paramValue = route.params.project_id
  const paramProjectId = Array.isArray(paramValue) ? paramValue[0] : paramValue

  const queryValue = route.query.project_id
  const queryProjectId = Array.isArray(queryValue) ? queryValue[0] : queryValue

  const raw = paramProjectId ?? queryProjectId
  const parsed = Number(raw)
  return Number.isInteger(parsed) && parsed > 0 ? parsed : null
})

const projectIdForQuery = computed<number>(() => projectId.value ?? 0)
const projectQuery = useProjectQuery(projectIdForQuery)
const projectName = computed(() => projectQuery.data.value?.name ?? 'Unknown project')
const layoutContextTitle = computed(() => `Add layout ${oldBarcode.value || '-'} to project ${projectName.value}`)
const layoutContextDescription = computed(
  () => `Enter the new barcode to copy control layout ${oldBarcode.value || '-'} into project ${projectName.value}.`,
)

watch(
  projectId,
  async (id) => {
    if (!id) return
    await layoutStore.fetchControlPlates()
  },
  { immediate: true },
)

const handleClick = (barcode: string) => {
  oldBarcode.value = barcode
  addLayoutDialog.value = true
}

const closeDialog = () => {
  addLayoutDialog.value = false
  newBarcode.value = ''
  oldBarcode.value = ''
}

const addLayout = async () => {
  const barcode = newBarcode.value.trim()
  if (!barcode) {
    toast.add({
      title: 'Please enter a barcode',
      color: 'warning',
      duration: 2500,
    })
    return
  }

  await layoutStore.addControlLayout({
    barcode_new: barcode,
    barcode_old: oldBarcode.value,
    project_id: projectId.value,
  })

  await Promise.all([projectQuery.refetch(), layoutStore.fetchControlPlates()])

  await queryClient.invalidateQueries({ queryKey: PROJECTS_QUERY_KEY })

  closeDialog()
  toast.add({
    title: 'Layout added',
    color: 'success',
    duration: 2500,
  })
}

const onWellSelected = () => {
  // Well details are not part of the control-layout selection page.
}
</script>

<template>
  <section class="select_layout_container">
    <p v-if="!projectId" class="px-2 text-sm text-slate-600">Missing project_id. Use `/layout/123`.</p>

    <p
      v-else-if="projectQuery.isPending.value || layoutStore.isLoadingControlPlates"
      class="px-2 text-sm text-slate-600"
    >
      Loading control layouts...
    </p>

    <p v-else-if="projectQuery.error.value" class="px-2 text-sm text-red-600">
      {{ projectQuery.error.value instanceof Error ? projectQuery.error.value.message : 'Failed to load project.' }}
    </p>

    <template v-else-if="projectQuery.data.value">
      <div class="project_context">
        Add layout to project: <span class="project_name">{{ projectName }}</span>
      </div>

      <section
        v-for="(plate, index) in layoutStore.controlPlates"
        :key="`control_plate_${plate.id}`"
        class="plate_container"
      >
        <div class="title">Plate {{ index + 1 }}: {{ plate.barcode }}</div>
        <button type="button" class="layout_link" @click="handleClick(plate.barcode)">
          Add layout {{ projectQuery.data.value.name }} >>
        </button>

        <div class="mt-4 overflow-auto">
          <DynamicPlate :plate="plate" @well-selected="onWellSelected" />
        </div>
      </section>
    </template>

    <WavesModalWrapper
      :open="addLayoutDialog"
      :title="layoutContextTitle"
      :description="layoutContextDescription"
      @update:open="addLayoutDialog = $event"
    >
      <template #body>
        <BaseField
          v-model="newBarcode"
          label="Enter barcode"
          placeholder="New barcode"
          input-class="w-full px-4 py-3 bg-white/70 backdrop-blur-sm outline-none ring-offset-0 focus:ring-2 focus:ring-lime-500 shadow rounded-full"
        />
      </template>

      <template #footer>
        <BaseButton
          label="Add"
          :on-click="addLayout"
          :loading="layoutStore.isAddingControlLayout"
          variant="primary"
          size="sm"
          width="auto"
        />
        <BaseButton
          label="Cancel"
          :on-click="closeDialog"
          variant="secondary"
          size="sm"
          width="auto"
          :disabled="layoutStore.isAddingControlLayout"
        />
      </template>
    </WavesModalWrapper>
  </section>
</template>

<style scoped>
.select_layout_container {
  width: 80%;
  max-width: 80vw;
}

.plate_container {
  padding: 20px;
}

.title {
  font-size: 1.5em;
  font-weight: bold;
  margin-bottom: 10px;
  color: #0f172a;
}

.layout_link {
  cursor: pointer;
  color: #0f766e;
  font-weight: 500;
  background: transparent;
  border: 0;
  padding: 0;
}

.layout_link:hover {
  text-decoration: underline;
}

.project_context {
  margin: 0.75rem 0.5rem 0.5rem;
  color: #0f172a;
  font-size: 1rem;
  font-weight: 500;
}

.project_name {
  color: #0f766e;
  font-weight: 700;
}
</style>
