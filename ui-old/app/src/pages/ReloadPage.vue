<template>
  <div class="fullscreen text-white text-center q-pa-md flex flex-center">
    <q-spinner-rings color="grey" size="2em" />
  </div>
</template>

<script setup lang="ts">
import {onMounted} from 'vue'

onMounted(() => {
  window.location.reload()
})
</script>
