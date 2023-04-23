<script setup lang="ts">
import {useManagementStore} from 'stores/management'
import {onMounted, ref} from 'vue'
import SelectedPaths from 'components/management/SelectedPaths.vue'
import MapCommand from 'components/management/MapCommand.vue'
import {useI18n} from 'vue-i18n'

const managementStore = useManagementStore()

const {t} = useI18n()

onMounted(async () => {
  await managementStore.initialize()
})

const tab = ref<string>('map')
</script>

<template>
  <q-page class="q-px-md">
    <h5 class="text-primary">{{ t('management.management') }}</h5>

    <SelectedPaths />

    <div class="q-gutter-y-md q-mt-xl" style="max-width: 800px">
      <q-card>
        <q-tabs
          v-model="tab"
          dense
          class="text-grey"
          active-color="primary"
          indicator-color="primary"
          align="justify"
          narrow-indicator>
          <q-tab name="map" :label="t('management.echo')"></q-tab>
          <q-tab name="import" :label="t('management.import')"></q-tab>
          <q-tab name="export" :label="t('management.export')"></q-tab>
        </q-tabs>

        <q-separator></q-separator>

        <q-tab-panels v-model="tab" animated>
          <q-tab-panel name="map">
            <MapCommand />
          </q-tab-panel>

          <q-tab-panel name="import">
            <div class="text-h6">Import</div>
            Gallia est omnis divisa in partes tres, quarum unam incolunt Belgae, aliam Aquitani, tertiam, qui
            ipsorum lingua Celtae, nostra Galli appellantur.
          </q-tab-panel>

          <q-tab-panel name="export">
            <div class="text-h6">Export</div>
            Ille mi par esse deo videtur, ille, si fas est, superare divos, qui sedens adversus identidem te
            spectat et audit.
          </q-tab-panel>
        </q-tab-panels>
      </q-card>
    </div>

    <!--    <div>-->
    <!--      <ul>-->
    <!--        <directory-tree-item :item="dataDirectory"></directory-tree-item>-->
    <!--      </ul>-->
    <!--    </div>-->
  </q-page>
</template>
