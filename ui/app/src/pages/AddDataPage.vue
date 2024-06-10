<script setup lang="ts">
import {useRoute} from 'vue-router'
import {onMounted, ref} from 'vue'
import {useProjectStore} from 'src/stores/project'
import {PlateInfo} from 'src/components/models'
import {useQuasar} from 'quasar'

const $q = useQuasar()

const route = useRoute()
const projectStore = useProjectStore()

onMounted(async () => {
  const exp_id = Number(route.params.experiment_id)
  await projectStore.getPrefilledPlateInfo(exp_id)
})

const downloadCSV = () => {
  const columnNames = [
    'measurement_label',
    'measurement_timestamp',
    'plate_barcode',
    'lib_plate_barcode',
    'replicate',
    'cell_type',
    'condition',
  ]

  const headerString = columnNames.join(',') + '\n'

  const csvContent = projectStore.prefilledPlateInfo
    .map((row: PlateInfo) => {
      // Map each column name to its corresponding value in the row object
      return columnNames
        .map(fieldName => {
          const value = row[fieldName as keyof PlateInfo]
          if (Array.isArray(value)) {
            return `"${value.join(',')}"`
          }
          return `"${value}"`
        })
        .join(',')
    })
    .join('\n')

  const csv = headerString + csvContent
  const blob = new Blob([csv], {type: 'text/csv'})
  const url = window.URL.createObjectURL(blob)
  const a = document.createElement('a')
  a.href = url
  a.download = 'experiment_data.csv'
  a.click()

  // Clean up by revoking the created URL
  window.URL.revokeObjectURL(url)
}
const savePlateInfo = async () => {
  alert('savePlateInfo')
  const result = await projectStore.savePlateInfo(
    Number(route.params.experiment_id),
    projectStore.prefilledPlateInfo
  )
  if (result === 'success') {
    $q.notify({
      message: 'Data saved successfully',
      color: 'positive',
      position: 'top',
      timeout: 2000,
    })
  } else {
    $q.notify({
      message: 'Error saving data',
      color: 'negative',
      position: 'top',
      timeout: 2000,
    })
  }
}
</script>

<template>
  <q-page class="q-px-md">
    <div class="text-h5 q-mt-lg q-mb-md q-pl-xl text-primary text-center">Add Experiment Data</div>

    <q-card class="my-card" flat>
      <q-card-section>
        <div class="text-body1">
          <vue-excel-editor
            v-model="projectStore.prefilledPlateInfo"
            :page="20"
            :no-paging="false"
            :filter-row="true"
            :autocomplete="true"
            :autocomplete-count="100"
            :readonly="false"
            :width="'100%'"
            :spellcheck="true"
            :new-if-bottom="true"
            :remember="true"
            :enter-to-south="true"
            :disable-panel-setting="false"
            :disable-panel-filter="false"
            :no-mouse-scroll="false">
            <vue-excel-column
              field="measurement_label"
              label="Measurement Label"
              type="string"
              width="150px" />
            <vue-excel-column
              field="measurement_timestamp"
              label="Measurement Timestamps"
              type="string"
              width="200px" />
            <vue-excel-column field="plate_barcode" label="Plate Barcode" type="string" width="120px" />
            <vue-excel-column
              field="lib_plate_barcode"
              label="Library Plate Barcode"
              type="string"
              width="150px" />
            <vue-excel-column field="replicate" label="Replicate" type="string" width="100px" />
            <vue-excel-column field="cell_type" label="Cell Type" type="string" width="100px" />
            <vue-excel-column field="condition" label="Condition" type="string" width="100px" />
          </vue-excel-editor>
        </div>
      </q-card-section>
      <q-card-actions align="left">
        <q-btn color="primary" label="Save" @click="savePlateInfo" />
        <q-btn color="primary" label="Download CSV" @click="downloadCSV"></q-btn>
      </q-card-actions>
    </q-card>
  </q-page>
</template>

<style scoped lang="sass">

.my-card
  max-width: 1000px
  margin: 0 auto
</style>
