<script setup lang="ts">
import {Project} from 'src/components/models'
import {defineProps, onMounted, PropType, ref} from 'vue'
import {useI18n} from 'vue-i18n'
import {api} from 'src/boot/axios'
import {useProjectStore} from 'stores/project'
import bus from 'src/eventBus'

const {t} = useI18n()
const projectStore = useProjectStore()

onMounted(() => {
  if (props.project) {
    experimentOptions.value = props.project.project.experiments.map(experiment => experiment.name)
  }
})

const props = defineProps({
  project: {
    type: Object as PropType<{project: Project}>,
    required: true,
  },
})

const selectedPlates = ref<string[]>([])
const selectedExperiment = ref<string | null>(null)
const experimentOptions = ref<string[]>([])

const movePlates = async () => {
  const payload = {
    plate_barcodes: selectedPlates.value,
    experiment: selectedExperiment.value,
  }

  try {
    const res = await api.post('/api/experiments/move_plates/', payload)
    console.log(res)
    await projectStore.initialize()
    bus.emit('project-updated')
  } catch (err) {
    console.error(err)
  }
}
</script>

<template>
  <section class="flex flex-row wrap">
    <article class="expansion-items-container">
      <q-expansion-item
        :label="exp.name"
        expand-separator
        v-for="exp in props.project.project.experiments"
        :key="`${exp.id}_exp`">
        <div v-for="plate in exp.plates" :key="`${plate.id}_plt`">
          <q-checkbox v-model="selectedPlates" :val="plate.barcode" :label="plate.barcode" />
        </div>
      </q-expansion-item>
    </article>
    <article class="q-ml-xl">
      <div>
        <q-select
          v-model="selectedExperiment"
          :options="experimentOptions"
          :label="t('label.chose_experiment')" />
        <q-btn @click="movePlates" class="q-mt-lg" color="secondary">
          {{ t('label.move_plates_to_experiment') }}
        </q-btn>
      </div>
    </article>
  </section>
</template>

<style scoped>
.expansion-items-container {
  width: 30%;
}
</style>
