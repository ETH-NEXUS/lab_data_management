<script setup lang="ts">
import {onMounted, ref} from 'vue'
import {api} from 'boot/axios'
import {useQuasar} from 'quasar'

const $q = useQuasar()
const thresholdUrl = '/api/thresholds/'
const redFlagUrl = '/api/compoundlib/redflag/'
const recalculateStatusUrl = '/api/compoundlib/recalculate_status/'

const dialog = ref(false)

onMounted(async () => {
  await getInfo()
  await getThresholds()
})
interface Threshold {
  id: number
  dmso: number
  amount: number
}

const thresholds = ref<Threshold>({dmso: 80, amount: 2.5, id: 1})
const redFlagInfo = ref<any>([])

const getInfo = async () => {
  try {
    const resp = await api.get(redFlagUrl)
    redFlagInfo.value = resp.data
  } catch (err) {
    console.error(err)
  }
}

const getThresholds = async () => {
  try {
    const resp = await api.get(thresholdUrl)
    thresholds.value = resp.data.results[0]
  } catch (err) {
    console.error(err)
  }
}

const updateThreshold = async () => {
  try {
    await api.patch(`${thresholdUrl}${thresholds.value.id}/`, {
      dmso: thresholds.value.dmso,
      amount: thresholds.value.amount,
    })
    dialog.value = false
    await getThresholds()
  } catch (err) {
    console.error(err)
  }
}

const recalculateStatus = async () => {
  try {
    $q.loading.show({
      message: 'Recalculating status of the library plates. Please wait, it can take some time...',
    })
    await api.get(recalculateStatusUrl)
    await getInfo()
    $q.loading.hide()
  } catch (err) {
    console.error(err)
    $q.loading.hide()
  }
}
</script>

<template>
  <q-page class="q-px-md">
    <div class="text-h5 q-mt-lg q-mb-md q-pl-xl text-primary text-center">Messages and Warnings</div>

    <q-card class="my-card" flat>
      <q-card-section class="q-pt-xs">
        <div class="text-body1 q-pl-md text-container">
          <div class="text-h6 q-mt-lg q-mb-md q-pl-sml text-primary">Empty Wells</div>
          <p class="text-body1 text-grey-8">
            Here, you can see a list of plates with problematic wells (i.e., wells that have an amount or DMSO
            concentration less than the threshold). This information is also represented on the navigation
            tree as a warning icon. You can change the threshold here and recalculate the status of the wells
            and plates.
          </p>
          <p class="text-body1 text-grey-8">
            Note that this information is derived from echo-transfer files; therefore, it only considers those
            plates that have been recently imported to the LDM via the echo-import process. For older
            projects, you can always resubmit the echo-import in the management tab to ensure that the
            calculation of the well status affects the plates.
          </p>
        </div>
      </q-card-section>
      <q-card-section>
        <div class="text-body1 q-pl-md text-container">
          <div class="text-h6 q-mb-md q-pl-sml text-primary">Problematic plates</div>
          <q-expansion-item
            v-for="(plates, libraryName) in redFlagInfo"
            :key="libraryName"
            :label="libraryName"
            class="text-primary exp">
            <q-expansion-item
              v-for="(wells, plateBarcode) in plates"
              :key="plateBarcode"
              :label="plateBarcode"
              class="text-secondary exp2">
              <q-list>
                <q-item v-for="well in wells" :key="well">
                  <q-item-section>
                    <span class="text-grey-7">{{ well }}</span>
                  </q-item-section>
                </q-item>
              </q-list>
            </q-expansion-item>
          </q-expansion-item>
        </div>
      </q-card-section>
      <q-card-section>
        <div class="text-body1 q-pl-md text-container">
          <div class="text-h6 q-mb-md q-pl-sml text-primary">Current Thresholds</div>
          <p class="text-body1 text-grey-8">
            <b>DMSO</b>
            : {{ thresholds.dmso }}%
          </p>
          <p class="text-body1 text-grey-8">
            <b>Amount</b>
            : {{ thresholds.amount }}µL
          </p>
          <q-btn label="Edit Threshold" color="primary" @click="dialog = true"></q-btn>
          <br />
          <q-btn
            label="Recalculate based on the current threshold"
            color="secondary"
            class="q-mt-md"
            @click="recalculateStatus"></q-btn>
        </div>
      </q-card-section>
    </q-card>
    <q-dialog v-model="dialog">
      <q-card>
        <q-card-section class="row q-pt-md q-pb-md">
          <div class="q-mb-md full-width">
            <q-input v-model="thresholds.dmso" label="DMSO (%)" type="number"></q-input>
            <q-input v-model="thresholds.amount" label="Amount (µL)" type="number" class="q-mt-md"></q-input>
          </div>
        </q-card-section>
        <q-card-actions align="right">
          <q-btn flat label="Submit" color="primary" @click="updateThreshold"></q-btn>
          <q-btn flat label="Cancel" color="negative" @click="dialog = false"></q-btn>
        </q-card-actions>
      </q-card>
    </q-dialog>
  </q-page>
</template>

<style scoped lang="sass">

.text-container
  max-width: 600px
  overflow-wrap: anywhere
  text-align: justify

.exp
  border: 1px solid #eee6e6
  border-radius: 2px
  background-color: #f0f0f0
  min-width: 400px

.exp2
  background-color: white
</style>
