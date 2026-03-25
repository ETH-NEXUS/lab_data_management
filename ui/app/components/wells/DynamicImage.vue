<script setup lang="ts">
import { onMounted, ref } from 'vue'
import { useAPI } from '~/composables/useAPI'

const props = withDefaults(
  defineProps<{
    url: string
    width?: string
    height?: string
  }>(),
  {
    width: '',
    height: '',
  },
)

const src = ref<string>()

onMounted(async () => {
  try {
    const endpoint = props.url.replace(/^\/api\//, '').replace(/^\//, '')
    const { data, error } = await useAPI<{ src?: string }>(endpoint, { method: 'GET' })
    if (!error.value) {
      src.value = data.value?.src
    }
  } catch (err) {
    console.error(err)
  }
})
</script>

<template>
  <img v-if="src" :src="src" :width="props.width" :height="props.height" alt="" />
  <span v-else class="inline-block h-8 w-8 animate-spin rounded-full border-2 border-teal-800/30 border-t-teal-800" />
</template>
