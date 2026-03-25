<script setup lang="ts">
import ldmLogo from '../assets/ldm_logo.png'
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
  width: 60%
  height: 60%
  display: flex
  align-items: center
  justify-content: center


.image
  max-width: 100%
  max-height: 100%
</style>
