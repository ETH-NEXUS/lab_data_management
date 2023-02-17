<script setup lang="ts">
import {useI18n} from 'vue-i18n'
import {defineProps, ref} from 'vue'
import {useQuasar} from 'quasar'
import {sidesData} from 'components/data'

import {useProjectStore} from 'stores/project'

const emit = defineEmits(['update'])

const props = defineProps({
  experimentId: {
    type: Number,
    required: true,
  },
})

const projectStore = useProjectStore()

const prefix = ref<string>('')
const number_of_plates = ref<number | null>(null)
const sides = ref<string[]>([])

const {t} = useI18n()
const $q = useQuasar()

const onCheck = () => {
  sides.value = []
  for (const side of sidesData) {
    const checkbox = document.getElementById(side.id) as HTMLInputElement
    if (checkbox.checked) {
      sides.value.push(side.label)
    }
  }
}

const onSubmit = async () => {
  if (sides.value.length === 0) {
    $q.notify({
      message: 'Please select at least one side',
      color: 'secondary',
    })
    return
  }

  if (prefix.value && number_of_plates.value && sides.value) {
    await projectStore.generateBarcodes(props.experimentId, prefix.value, number_of_plates.value, sides.value)
    emit('update')
  }
}

const onReset = () => {
  prefix.value = ''
  number_of_plates.value = 0
}
</script>

<template>
  <q-card class="formDialog">
    <q-card-section class="row items-center q-pb-none">
      <div class="text-h6">{{ t('action.generate_barcodes') }}</div>
      <!--      <q-space></q-space>-->
      <!--      <q-btn icon="close" flat round dense v-close-popup></q-btn>-->
    </q-card-section>

    <q-card-section>
      <div class="q-pa-md" style="max-width: 400px">
        <q-form @submit="onSubmit" @reset="onReset" class="q-gutter-md full-width">
          <q-input
            v-model="prefix"
            :label="t('action.enter_prefix')"
            lazy-rules
            :rules="[val => (val && val.length > 0) || t('action.validation_prefix')]"></q-input>

          <q-input
            type="number"
            v-model="number_of_plates"
            :label="t('action.enter_number_of_plates')"
            lazy-rules
            :rules="[
              val => (val !== null && val !== '') || t('action.validation_number_of_plates'),
              val => val > 0 || t('action.validation_positive_number'),
            ]"></q-input>

          <div class="outerDiv">
            <div class="innerDiv">
              <div v-for="side in sidesData" :key="side.id">
                <input
                  type="checkbox"
                  :id="side.id"
                  :style="side.styleCheckbox"
                  :value="side.id"
                  @change="onCheck()" />
                <label :for="side.id" :style="side.styleLabel">{{ side.label }}</label>
              </div>
            </div>
          </div>
          <q-separator class="q-mt-lg"></q-separator>
          <div class="q-mt-lg">
            <q-btn label="Submit" type="submit" color="primary"></q-btn>
            <q-btn label="Reset" type="reset" color="primary" flat class="q-ml-sm"></q-btn>
            <q-btn flat label="Cancel" v-close-popup></q-btn>
          </div>
        </q-form>
      </div>
    </q-card-section>
  </q-card>
</template>

<style>
.formDialog {
  width: 35vmin;
}

.outerDiv {
  margin-top: 20px;
  padding: 20px 20px 50px 75px;
}
.innerDiv {
  width: 100px;
  height: 60px;
  background-color: #4ca4f3;
  position: relative;
  border-radius: 5px;
}
</style>
