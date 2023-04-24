<script setup lang="ts">
import DynamicForm from './DynamicForm.vue'
import {Options, FormData} from '../models'
import {useManagementStore} from 'stores/management'
import {onMounted} from 'vue'
import {useQuasar} from 'quasar'
import {useI18n} from 'vue-i18n'

const props = defineProps<{
  options: Options
  command: string
}>()

const managementStore = useManagementStore()
const $q = useQuasar()
const {t} = useI18n()

onMounted(() => {
  managementStore.commandOutput = ''
})

const onSubmit = async (formData: FormData) => {
  formData['command'] = props.command
  // $q.loading.show({
  //   message: t('info.running_in_progress'),
  // })
  await managementStore.runCommand(formData)
  //$q.loading.hide()
}
</script>

<template>
  <div>
    <dynamic-form :options="options" @submit="onSubmit" />
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
