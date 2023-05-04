<script setup lang="ts">
import {useI18n} from 'vue-i18n'
import ldmLogo from '../assets/ldm_logo.png'

const {t} = useI18n()
</script>

<template>
  <q-page class="row items-center justify-evenly">
    <div class="image-container">
      <img class="image" :src="ldmLogo" alt="LDM Logo" />
    </div>
  </q-page>
</template>

<style scoped lang="sass">
.image-container
  width: 80%
  height: 80%
  display: flex
  align-items: center
  justify-content: center


.image
  max-width: 100%
  max-height: 100%
</style>
