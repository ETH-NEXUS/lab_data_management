<script setup lang="ts">
import { ref, defineProps, onMounted } from 'vue'
import { api} from '../boot/axios'

const props = defineProps({
  url: {
    type: String,
    required: true
  },
  width: {
    type: String,
    default: ''
  },
  height: {
    type: String,
    default: ''
  }
})

const src = ref<string>()

onMounted(async () => {
  try {
    const resp = await api.get(props.url)
    src.value = resp.data.src
  } catch (err) {
    console.error(err)
  }
})
</script>

<template>
  <q-img v-if="src" :src="src" :width="props.width" :height="props.height"/>
  <q-spinner-grid v-else color="primary" size="2em" />
</template>
