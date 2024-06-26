<script setup lang="ts">
import {onMounted, ref} from 'vue'
import {useProjectStore} from 'stores/project'
import DynamicPlate from 'components/plate/DynamicPlate.vue'
import {Project} from 'components/models'
import {useRoute} from 'vue-router'
import {api} from 'src/boot/axios'
import {useI18n} from 'vue-i18n'
import bus from 'src/eventBus'
import {useQuasar} from 'quasar'

const $q = useQuasar()
const {t} = useI18n()

const route = useRoute()
const projectStore = useProjectStore()
const project = ref<Project | null>(null)
const addLayoutDialog = ref<boolean>(false)
const newBarcode = ref<string>('')
const oldBarcode = ref<string>('')

onMounted(async () => {
  if (projectStore.projects.length === 0) {
    await projectStore.initialize()
  }

  project.value = await getProject()
  await projectStore.getControlPlates()
})

const getProject = async () => {
  const res = await api.get(`/api/projects/${route.params.project_id}`)
  return res.data
}

const addLayout = async () => {
  if (!newBarcode.value) {
    $q.notify({
      message: 'Please enter a barcode',
      color: 'warning',
      position: 'top',
    })
    return
  }
  const payload = {
    barcode_new: newBarcode.value,
    barcode_old: oldBarcode.value,
    project_id: project.value?.id,
  }

  const res = await api.post('/api/add_control_layout/', payload, {
    headers: {
      'Content-Type': 'application/json',
    },
  })
  console.log(res.status)
  bus.emit('project-updated')
}

const handleClick = (barcode: string) => {
  oldBarcode.value = barcode
  addLayoutDialog.value = true
}
</script>

<template>
  <div class="select_layout_container" v-if="project">
    <section
      class="plate_container"
      v-for="(plate, index) in projectStore.controlPlates"
      :key="`control_plate_${plate.id}`">
      <div class="text-primary title">Plate {{ index + 1 }}: {{ plate.barcode }}</div>
      <div @click="handleClick(plate.barcode)" class="cursor-pointer text-secondary">
        {{ t('action.add_layout') }} {{ project.name }} >>
      </div>
      <DynamicPlate :plate="plate" :min="0" :max="1" />
    </section>
    <q-dialog v-model="addLayoutDialog">
      <q-card class="dialog-card">
        <q-card-section>
          <q-input v-model="newBarcode" label="Enter barcode"></q-input>
        </q-card-section>
        <q-card-actions align="right">
          <q-btn flat label="Add" color="primary" @click="addLayout" />
          <q-btn flat label="Cancel" color="primary" @click="addLayoutDialog = false" />
        </q-card-actions>
      </q-card>
    </q-dialog>
  </div>
</template>

<style scoped lang="sass">

.select_layout_container
  width: 80%
  max-width: 80vw

.plate_container
  padding: 20px

.title
    font-size: 1.5em
    font-weight: bold
    margin-bottom: 10px

.dialog-card
    width: 400px
</style>
