<script setup lang="ts">
import {api} from 'src/boot/axios'
import {ref, onMounted} from 'vue'
import {useProjectStore} from 'stores/project'
import {useQuasar} from 'quasar'
import {useI18n} from 'vue-i18n'

const props = defineProps({
  experimentName: {
    type: String,
    required: true,
  },
  label: {
    type: String,
    required: true,
  },
})

onMounted(async () => {
  await getInputNotebookOptions()
})

const notebookOptions = ref<string[]>([])
const projectStore = useProjectStore()
const $q = useQuasar()
const {t} = useI18n()
const inputNotebookPath = ref<string | null>(null)

const getInputNotebookOptions = async () => {
  const response = await api.post('/api/list_files/', {
    notebooks_dir: '/notebooks/input',
    file_format: '.ipynb',
  })
  console.log('input', response.data)
  notebookOptions.value = response.data.notebooks
}

const submit = async () => {
  if (!inputNotebookPath.value) {
    $q.notify({
      message: 'Please select a notebook template  to generate a report from.',
      type: 'negative',
    })
  } else {
    $q.loading.show({
      message: t('info.generation_in_progress'),
    })
    await projectStore.generateReport(props.experimentName, props.label)
    $q.loading.hide()
  }
}
</script>

<template>
  <q-card class="card">
    <q-card-section>
      <q-select v-model="inputNotebookPath" :options="notebookOptions" label="Notebook Template">
        <template v-slot:option="scope">
          <q-item v-bind="scope.itemProps" v-on="scope.itemEvents">
            <q-item-section>
              {{ scope.opt.split('/').pop() }}
            </q-item-section>
          </q-item>
        </template>
      </q-select>
    </q-card-section>
    <q-card-actions align="right">
      <q-btn flat label="Cancel" color="primary" v-close-popup />
      <q-btn flat label="Generate" color="primary" v-close-popup @click="submit" />
    </q-card-actions>
  </q-card>
</template>

<style scoped lang="sass">
.card
  width: 400px
  max-width: 90vw
  min-height: 150px
  max-height: 90vh
  overflow: auto
</style>
