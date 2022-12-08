<script setup lang="ts">
import {Plate} from 'src/components/models'
import {ref, onMounted} from 'vue'
import {api} from '../boot/axios'
import {useRoute} from 'vue-router'
import DynamicPlate from '../components/DynamicPlate.vue'
import {handleError} from '../helpers/errorHandling'
import {Well} from '../components/models'
import WellInfo from '../components/WellInfo.vue'

const route = useRoute()

const plate = ref<Plate | null>(null)
const splitter = ref<number>(50)
const selectedWell = ref<Well>()

onMounted(async () => {
    try {
        const resp = await api.get(`/api/plates/?barcode=${route.params.barcode}`)
        if (resp.data.results.length === 1) {
            plate.value = resp.data.results[0]
        } else if (resp.data.results.length === 0) {
            handleError(`No plate found with barcode ${route.params.barcode}.`)
        } else {
            handleError(`Multiple plates found with barcode ${route.params.barcode}.`)
        }
    } catch (err) {
        handleError(err)
    }
})
</script>

<template>
    <q-page class="row items-top q-px-md" :key="`${route.params.barcode}`">
        <q-splitter v-model="splitter" class="full-width">
            <template v-slot:before>
                <div v-if="plate">
                    <h2>{{ plate.barcode }}</h2>
                    <dynamic-plate :plate="plate" @well-selected="well => (selectedWell = well)" />
                </div>
                <q-spinner-grid v-else color="primary" size="5em" class="absolute-center" />
            </template>

            <template v-slot:after>
                <div v-if="selectedWell" class="q-px-md">
                    <h2>{{ selectedWell.hr_position }}</h2>
                    <well-info :well="selectedWell" />
                </div>
            </template>
        </q-splitter>
    </q-page>
</template>

<style scoped lang="sass">
h2
  font-family: 'Courier New', Courier, monospace
</style>
