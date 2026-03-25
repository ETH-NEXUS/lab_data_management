<script setup lang="ts">
import DynamicForm from './DynamicForm.vue'
import {Options, GeneralFormData} from '../models'
import {useManagementStore} from 'stores/management'
import {onMounted} from 'vue'
import {useUserStore} from 'stores/user'
import {storeToRefs} from 'pinia'
import bus from 'src/eventBus'

const props = defineProps<{
  options: Options
  command: string
  what: string
}>()

const managementStore = useManagementStore()

const {user} = storeToRefs(useUserStore())

onMounted(() => {
  managementStore.commandOutput = ''
})

const onSubmit = async (formData: GeneralFormData) => {
  managementStore.commandOutput = `Executing command: ${props.command}\n...`
  formData['room_name'] = user.value?.id.toString() || 'room_name'
  formData['command'] = props.command
  if (props.what && props.command === 'import') {
    if (props.what === 'control_plate') {
      formData['what'] = 'library_plate'
      formData['is_control_plate'] = true
    } else if (props.what === 'library_plate') {
      formData['what'] = 'library_plate'
      formData['is_control_plate'] = false
    } else {
      formData['what'] = props.what
    }
  }

  await managementStore.runCommand(formData)
  bus.emit('management-command')
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
