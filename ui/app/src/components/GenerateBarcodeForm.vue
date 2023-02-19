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
  prefilledData: {
    type: Object,
    required: false,
  },
  edit: {
    type: Boolean,
    required: false,
    default: false,
  },
})

const projectStore = useProjectStore()

const prefix = ref<string>(props.prefilledData?.prefix || '')
const number_of_plates = ref<number | null>(props.prefilledData?.number_of_plates || null)
const sides = ref<string[]>(props.prefilledData?.sides || [])

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
    if (props.prefilledData?.id) {
      await projectStore.updateBarcode(
        props.prefilledData.id,
        prefix.value,
        number_of_plates.value,
        sides.value
      )
      emit('update')
    } else {
      await projectStore.generateBarcodes(
        props.experimentId,
        prefix.value,
        number_of_plates.value,
        sides.value
      )
      emit('update')
      return
    }
  }
}

const onReset = () => {
  prefix.value = ''
  number_of_plates.value = 0
}
</script>

<template>
  <q-card :class="`${props.edit ? 'editFormDialog' : 'formDialog'} q-mb-lg`" :bordered="false">
    <q-card-section class="row items-center q-pb-none">
      <div class="text-h6">
        {{ props.edit ? 'Edit barcode specifications' : t('action.generate_barcodes') }}
      </div>
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

          <div class="outerDiv" v-if="!props.edit">
            <div class="innerDiv">
              <div v-for="side in sidesData" :key="side.id">
                <input
                  type="checkbox"
                  :id="side.id"
                  :style="side.styleCheckbox"
                  :value="side.id"
                  :checked="props.prefilledData?.sides.includes(side.label)"
                  @change="onCheck()" />
                <label :for="side.id" :style="side.styleLabel">
                  {{ side.label }}
                </label>
              </div>
            </div>
          </div>
          <div v-else>
            <q-select
              v-model="sides"
              :label="t('action.select_sides')"
              :options="sidesData.map(side => side.label)"
              multiple
              use-chips
              stack-label
              :rules="[val => (val && val.length > 0) || t('action.validation_sides')]"
              :emit-value="true"
              :map-options="true"
              :fill-input="true"
              :input-debounce="0"
              :max-height="200"
              :style="{width: '100%'}"></q-select>
          </div>

          <div class="q-mt-lg">
            <q-btn label="Submit" type="submit" color="primary"></q-btn>
            <q-btn label="Reset" type="reset" color="primary" flat class="q-ml-sm"></q-btn>
          </div>
        </q-form>
      </div>
    </q-card-section>
  </q-card>
</template>

<style>
.formDialog {
  width: 20%;
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
