<script setup lang="ts">
import {defineProps} from 'vue'
import {useI18n} from 'vue-i18n'
import {WellDetails} from 'components/models'
import {Well} from 'components/models'

const props = defineProps({
  well: {
    type: Object as () => Well | WellDetails | undefined,
    required: true,
  },
  row: {
    type: Number,
    required: true,
  },
  col: {
    type: Number,
    required: true,
  },
  selectedMeasurement: {
    default: null,
  },
  selectedTimestampIdx: {
    type: Number,
    required: true,
  },
})

const {t} = useI18n()

const measurement = (well: WellDetails) => {
  if (props.selectedMeasurement && well.measurements) {
    if (props.selectedMeasurement in well.measurements) {
      if (well.measurements[props.selectedMeasurement].length > props.selectedTimestampIdx) {
        return well.measurements[props.selectedMeasurement][props.selectedTimestampIdx]
      }
    }
  }
  return null
}
</script>

<template>
  <q-tooltip class="tooltip" v-if="props.well" anchor="top middle" self="bottom middle" :offset="[5, 5]">
    <b>{{ well.hr_position }}</b>
    ({{ well.type }})
    <small>({{ well.position }})</small>
    <hr />
    <b v-if="well.compounds" class="q-mt-md">{{ t('label.compounds') }}</b>
    <table v-if="well.compounds">
      <tr v-for="compound in well.compounds" :key="compound">
        <b>{{ compound }}</b>
      </tr>
    </table>
    <b v-if="props.selectedMeasurement" class="q-mt-md">{{ t('label.measurements') }}</b>
    <br />
    <span v-if="props.selectedMeasurement" class="q-mt-md">
      {{ props.selectedMeasurement }}: {{ measurement(well) }}
    </span>
  </q-tooltip>
</template>

<style scoped></style>
