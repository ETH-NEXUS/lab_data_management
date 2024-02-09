<script setup lang="ts">
import {ref} from 'vue'
import {useProjectStore} from 'stores/project'
import {useQuasar} from 'quasar'
import {useI18n} from 'vue-i18n'

const props = defineProps({
  experimentName: {
    type: String,
    required: true,
  },

  labels: {
    type: Array,
    required: true,
  },
})

const projectStore = useProjectStore()
const $q = useQuasar()
const {t} = useI18n()
const selectedLabel = ref<string | null>(null)

const submit = async () => {
  if (!selectedLabel.value) {
    $q.notify({
      message: 'Please select a measurement label to download the data for.',
      type: 'negative',
    })
  } else {
    $q.loading.show({
      message: t('info.generation_in_progress'),
    })

    await projectStore.downloadCSVData(props.experimentName, selectedLabel.value)
    $q.loading.hide()
  }
}
</script>

<template>
  <q-card class="card">
    <q-card-section>
      <q-select
        class="q-mt-lg"
        label="Measurement label"
        :options="props.labels"
        v-model="selectedLabel"></q-select>
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
