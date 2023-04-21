<script setup lang="ts">
import {storeToRefs} from 'pinia'
import {useSettingsStore} from 'stores/settings'
import {useI18n} from 'vue-i18n'

const {platePage} = storeToRefs(useSettingsStore())
const {t} = useI18n()

const props = defineProps({
  showSquareCompoundType: {
    type: Boolean,
    default: true,
  },
})
</script>

<template>
  <q-toggle
    class="q-mb-md"
    size="sm"
    checked-icon="check"
    v-model="platePage.smallerMapView"
    :label="t('label.smaller_map_view')"
    right-label
    color="secondary"></q-toggle>
  <q-toggle
    v-if="showSquareCompoundType"
    class="q-mb-md"
    size="sm"
    checked-icon="check"
    v-model="platePage.squareCompoundType"
    :label="t('label.show_type_as_square')"
    right-label
    color="secondary"></q-toggle>
</template>

<style scoped></style>
