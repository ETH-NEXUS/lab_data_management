<script setup lang="ts">
import DynamicForm from './DynamicForm.vue'
import {Options, GeneralFormData} from '../models'
import {useManagementStore} from 'stores/management'
import {onMounted} from 'vue'
import {useQuasar} from 'quasar'
import {useI18n} from 'vue-i18n'
import {useUserStore} from 'stores/user'
import {storeToRefs} from 'pinia'

const props = defineProps<{
  options: Options
  command: string
}>()

const managementStore = useManagementStore()
const $q = useQuasar()
const {t} = useI18n()
const {user} = storeToRefs(useUserStore())

onMounted(() => {
  managementStore.commandOutput = ''
})

const onSubmit = async (formData: GeneralFormData) => {
  formData['room_name'] = user.value?.id.toString() || 'room_name'
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
      <div v-if="managementStore.commandOutput" class="terminal-window q-mt-lg">
        <pre>{{ managementStore.commandOutput }}</pre>
      </div>
    </div>
  </div>
</template>

<style scoped lang="sass">


.terminal-window
  background-color: #1e1e1e
  color: #f8f8f2
  font-family: 'Courier New', Courier, monospace
  font-size: 14px
  padding: 10px
  border-radius: 5px
  overflow: auto
  height: 200px
  max-width: 100%

  pre
    margin: 0
    padding: 0
    font-family: 'Courier New', Courier, monospace
    font-size: 14px
    white-space: pre-wrap
</style>
