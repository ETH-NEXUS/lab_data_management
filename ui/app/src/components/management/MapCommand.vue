<script setup lang="ts">
import DynamicForm from './DynamicForm.vue'
import {Options, FormData} from '../models'
import {useManagementStore} from 'stores/management'
import {ref, onMounted} from 'vue'
import {useQuasar} from 'quasar'
import {useI18n} from 'vue-i18n'

const managementStore = useManagementStore()
const $q = useQuasar()
const {t} = useI18n()

onMounted(() => {
  managementStore.commandOutput = ''
})

const mapCommandOptions: Options = {
  // probably there is a way to get this data dynamically from the backend
  machine: {
    type: 'str',
    choices: ['echo', 'm1000'],
    label: 'Machine to map from',
    required: true,
  },
  path: {
    type: 'str',
    label: 'Path to the directory containing the mapping files',
    required: true,
  },
  mapping_file: {
    type: 'str',
    label: 'For echo mapping only: a yml file with the column headers, otherwise default headers are used',
    required: false,
  },

  experiment_name: {
    type: 'str',
    label: 'Experiment name',
    required: true,
  },
  debug: {
    type: 'bool',
    label: 'Enable debug mode',
    required: false,
  },
}

const onSubmit = async (formData: FormData) => {
  formData['command'] = 'map'
  $q.loading.show({
    message: t('info.running_in_progress'),
  })
  await managementStore.runCommand(formData)
  $q.loading.hide()
}
</script>

<template>
  <div>
    <dynamic-form :options="mapCommandOptions" @submit="onSubmit" />
  </div>
  <div>
    <div>
      <div class="text-caption" v-if="managementStore.commandOutput">Logs:</div>
      <div v-if="managementStore.commandOutput">
        <pre>{{ managementStore.commandOutput }}</pre>
      </div>
    </div>
  </div>
</template>

<style scoped></style>
