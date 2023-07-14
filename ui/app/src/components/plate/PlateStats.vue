<script setup lang="ts">
import {defineProps} from 'vue'
import {useI18n} from 'vue-i18n'

const {t} = useI18n()

const props = defineProps({
  z_prime: {
    default: null,
  },
  ssmd: {
    default: null,
  },
})
</script>

<template>
  <div v-if="props.z_prime" class="q-mt-sm">
    <b>{{ t('label.z_prime') }}:</b>
    <span
      :class="{
        'measurement': true,
        'q-mx-xs': true,
        'q-pa-xs': true,
        'bg-green-3': z_prime >= 0.5 && z_prime <= 1,
        'bg-orange-3': z_prime >= 0 && z_prime < 0.5,
        'bg-red-3': z_prime < 0,
      }">
      {{ z_prime.toFixed(2) }}
    </span>
  </div>
  <div v-if="props.ssmd" class="q-mt-sm">
    <b>{{ t('label.ssmd') }}:</b>
    <span
      :class="{
        'measurement': true,
        'q-mx-xs': true,
        'q-pa-xs': true,
        'bg-green-3': props.ssmd >= 12,
        'bg-orange-3': props.ssmd >= 6 && props.ssmd < 12,
        'bg-red-3': props.ssmd < 6,
      }">
      {{ ssmd.toFixed(2) }}
    </span>
  </div>
</template>

<style scoped lang="sass"></style>
