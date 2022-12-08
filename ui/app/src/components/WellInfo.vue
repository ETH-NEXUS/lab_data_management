<script setup lang="ts">
import {defineProps, PropType} from 'vue'
import {Well} from './models'
import DynamicImage from './DynamicImage.vue'
import {useI18n} from 'vue-i18n'

const {t} = useI18n()

const props = defineProps({
    well: {
        type: Object as PropType<Well>,
        required: true,
    },
})
</script>

<template>
    <div class="container full-width">
        <div class="row">
            <div class="col-12">
                <h4>{{ t('title.compound') }}</h4>
            </div>
            <div class="col-12">
                <table>
                    <tr>
                        <th>{{ t('label.identifier') }}</th>
                        <td>{{ props.well.compound.identifier }}</td>
                    </tr>
                    <tr>
                        <th>{{ t('label.smile') }}</th>
                        <td>{{ props.well.compound.smile }}</td>
                    </tr>
                    <tr>
                        <td colspan="2">
                            <dynamic-image
                                :url="`/api/wells/${props.well.id}/structure/`"
                                width="250px"
                                :key="props.well.id" />
                        </td>
                    </tr>
                </table>
                <hr />
            </div>
            <div class="col-12">
                <h4>{{ t('title.well') }}</h4>
            </div>
            <div class="col-12">
                <table>
                    <tr>
                        <th>{{ t('label.amount') }}</th>
                        <td>{{ props.well.amount }}</td>
                    </tr>
                </table>
                <hr />
            </div>
            <div class="col-12">
                <h4>{{ t('title.measurements') }}</h4>
            </div>
            <div class="col-12">
                <table>
                    <thead>
                        <tr>
                            <th>{{ t('label.name') }}</th>
                            <th>{{ t('label.abbrev') }}</th>
                            <th>{{ t('label.value') }}</th>
                            <th>{{ t('label.unit') }}</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr v-for="measurement in props.well.measurements" :key="measurement.name">
                            <td>{{ measurement.name }}</td>
                            <td>{{ measurement.abbrev }}</td>
                            <td>{{ measurement.value }}</td>
                            <td>{{ measurement.unit }}</td>
                        </tr>
                    </tbody>
                </table>
                <hr />
            </div>
        </div>
    </div>
</template>

<style scoped lang="sass">
td, th
 text-align: left
h4
  font-size: 1.5em
  font-weight: bold
</style>
