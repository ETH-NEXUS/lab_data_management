<script setup lang="ts">
import {reactive, onMounted} from 'vue'
import {Options, GeneralFormData} from 'components/models'
import {useQuasar} from 'quasar'
import {useI18n} from 'vue-i18n'

const props = defineProps<{
  options: Options
}>()

const emit = defineEmits(['submit'])
const $q = useQuasar()

const {t} = useI18n()

const form = reactive({} as GeneralFormData)

onMounted(() => {
  initializeForm()
})

const initializeForm = () => {
  for (const key in props.options) {
    if (props.options[key].type === 'bool') {
      form[key] = false
    } else {
      form[key] = ''
    }
  }
}

const onSubmit = () => {
  if (Object.keys(form).includes('experiment_name') && form.experiment_name === '') {
    $q.notify({
      message: 'Please enter an experiment name',
      color: 'negative',
      position: 'top',
    })
  }
  if (Object.keys(form).includes('plate_barcode') && form.plate_barcode === '') {
    $q.notify({
      message: 'Please enter the barcode of the plate',
      color: 'negative',
      position: 'top',
    })
  }
  if (Object.keys(form).includes('library_name') && form.library_name === '') {
    $q.notify({
      message: 'Please enter the library name',
      color: 'negative',
      position: 'top',
    })
  }
  if (Object.keys(form).includes('project_name') && form.project_name === '') {
    $q.notify({
      message: 'Please enter the project name',
      color: 'negative',
      position: 'top',
    })
  }
  console.log('Form data:', form)
  emit('submit', form)
}
</script>

<template>
  <div class="q-pa-md formContainer">
    <q-form @submit="onSubmit" class="q-gutter-md">
      <div v-for="(option, key) in props.options" :key="option.label">
        <q-select
          v-if="option.type === 'str' && option.choices"
          v-model="form[key]"
          :options="option.choices"
          :hint="option.label"
          :required="option.required" />

        <q-input
          v-else-if="option.type === 'str'"
          v-model="form[key]"
          :hint="option.label"
          type="text"
          :required="option.required" />

        <q-toggle v-else-if="option.type === 'bool'" v-model="form[key]" :label="option.label" />
      </div>

      <q-btn color="secondary" type="submit" :label="t('action.submit')" class="q-mt-lg" />
    </q-form>
  </div>
</template>

<style>
.formContainer {
  max-width: 700px;
}
</style>
