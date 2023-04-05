<script setup lang="ts">
import {ref, defineProps, onMounted, defineEmits} from 'vue'
import {useQuasar} from 'quasar'

const props = defineProps({
  labels: {
    type: Array<string>,
    required: true,
  },
})

const emit = defineEmits(['calculate'])

const $q = useQuasar()

onMounted(async () => {
  console.log('mounted')
})

const buttons: string[] = ['1', '2', '3', '+', '4', '5', '6', '-', '7', '8', '9', '*', 'C', '0', '/', 'ln']
const currentExpression = ref<string>('')
const previousButton = ref<string>('')
const newLabel = ref<string>('')
const usedLabels = ref<string[]>([])

const handleButtonClick = (button: string) => {
  switch (button) {
    case 'C':
      currentExpression.value = ''
      break
    case 'ln':
      currentExpression.value += 'ln('
      previousButton.value = button
      break
    default:
      if (previousButton.value === 'ln') {
        currentExpression.value += `${button})`
      } else {
        currentExpression.value += button
      }
      previousButton.value = button
      break
  }
}
const handleDragStart = (label: string, event: DragEvent) => {
  event.dataTransfer?.setData('text/plain', label)
}
const handleDrop = (event: DragEvent) => {
  const label = event.dataTransfer?.getData('text/plain')
  currentExpression.value += label
  if (previousButton.value === 'ln') {
    currentExpression.value += ')'
  }
  if (label) {
    usedLabels.value.push(label)
  }
}

const addMeasurement = () => {
  if (!newLabel.value) {
    $q.notify({
      type: 'negative',
      message: 'Please enter a name for the new measurement',
    })

    return
  }
  if (props.labels?.includes(newLabel.value)) {
    $q.notify({
      type: 'negative',
      message: 'Please use a different name for the new measurement',
    })

    return
  }

  // check if used labels contain any of the values in labels
  const usedLabelsSet = new Set(usedLabels.value)
  const labelsSet = new Set(props.labels)
  const intersection = new Set([...usedLabelsSet].filter(x => labelsSet.has(x)))

  // the measurement will not be added if no labels were used
  if (intersection.size === 0) {
    $q.notify({
      type: 'negative',
      message: 'Please use at least one existing measurement',
    })

    return
  }

  emit('calculate', currentExpression.value, newLabel.value, usedLabels.value)
}
</script>

<template>
  <div>
    <div class="calculator">
      <input
        type="text"
        v-model="currentExpression"
        class="calculator__input"
        @drop.prevent="handleDrop($event)"
        @dragover.prevent />
      <div class="calculator__buttons">
        <button
          v-for="button in buttons"
          :key="button"
          @click="handleButtonClick(button)"
          class="calculator__button">
          {{ button }}
        </button>
      </div>
    </div>
    <div class="labels">
      <div
        v-for="(label, index) in labels"
        :key="index"
        class="label"
        :draggable="true"
        @dragstart="handleDragStart(label, $event)">
        {{ label }}
      </div>
    </div>

    <div class="q-pa-md q-gutter-md">
      <div class="col-3">
        <q-input v-model="newLabel" type="text" hint="New measurement name" outlined></q-input>
      </div>
      <div class="col-6 text-right">
        <q-btn icon="calculate" color="secondary" class="q-ml-xs" @click="addMeasurement">Apply</q-btn>
      </div>
    </div>
  </div>
</template>

<style lang="sass">


.calculator
  display: flex
  flex-direction: column
  align-items: center
  padding: 20px
  margin: 50px
  background-color: #f2f2f2
  border-radius: 10px

.calculator__input
  margin-bottom: 10px
  padding: 10px
  font-size: 14px
  text-align: right
  width: 100%
  box-sizing: border-box

.calculator__buttons
  display: grid
  grid-template-columns: repeat(4, 1fr)
  grid-gap: 10px

.calculator__button
  background-color: #ffffff
  color: #000000
  border: none
  border-radius: 5px
  padding: 10px
  font-size: 18px
  cursor: pointer
  transition: all 0.2s ease-in-out

.calculator__button:hover
  background-color: #cccccc

.labels
  display: flex
  flex-wrap: wrap
  justify-content: center
  align-items: center
  margin-top: 20px

.label
  padding: 10px 20px
  margin: 10px
  border-radius: 20px
  background-color: #f2f2f2
  font-size: 12px
  font-weight: bold
  color: #333
  box-shadow: 2px 2px 4px rgba(0, 0, 0, 0.2)
  cursor: move

.label:hover
  background-color: #cccccc

.label:active
  box-shadow: none
  background-color: rgba(255, 255, 255, 0.5)
</style>
