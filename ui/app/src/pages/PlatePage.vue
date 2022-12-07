<script setup lang="ts">
import {Plate} from 'src/components/models'
import {ref, onMounted} from 'vue'
import {api} from '../boot/axios'
import {useRoute} from 'vue-router'
import DynamicPlate from '../components/DynamicPlate.vue'
import {handleError} from '../helpers/errorHandling'

const route = useRoute()

const plate = ref<Plate | null>(null)

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
    <q-page class="row items-center justify-evenly" :key="`${route.params.barcode}`">
        <div v-if="plate">
            <h2>{{ plate.barcode }}</h2>
            <dynamic-plate :plate="plate" />
        </div>
        <q-spinner-grid v-else color="primary" size="5em" />
    </q-page>
</template>
