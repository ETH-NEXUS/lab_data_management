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
    <q-card-section>
      <div style="max-width: 400px">
        <div class="text-h6 q-pb-md">
          {{ props.edit ? 'Edit barcode specifications' : t('action.generate_barcodes') }}
        </div>
        <q-form dense @submit="onSubmit" @reset="onReset" class="full-width">
          <p class="text-body2">{{ t('action.enter_prefix') }}:</p>
          <q-input
            :dense="true"
            v-model="prefix"
            lazy-rules
            :rules="[val => (val && val.length > 0) || t('action.validation_prefix')]"></q-input>

          <p class="text-body2">{{ t('action.enter_number_of_plates') }}:</p>

          <q-input
            :dense="true"
            type="number"
            v-model="number_of_plates"
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

          <div class="q-mt-lg fir row wrap justify-end">
            <q-btn flat label="Reset" type="reset" color="primary"></q-btn>
            <q-btn class="q-ml-sm" flat label="Ok" type="submit" color="primary"></q-btn>
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
  padding: 20px 20px 50px 60px;
}
.innerDiv {
  width: 100px;
  height: 60px;
  background-color: #4ca4f3;
  position: relative;
  border-radius: 5px;
}
</style>
