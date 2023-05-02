<script setup lang="ts">
import {useI18n} from 'vue-i18n'

const {t} = useI18n()
</script>

<template>
  <q-page class="row items-center justify-evenly">
    <h1>{{ t('title.app') }}</h1>
  </q-page>
</template>
