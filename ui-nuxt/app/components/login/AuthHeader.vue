<script setup lang="ts">
type Props = {
  logoSrc: string
  subtitle: string
}

const props = defineProps<Props>()
</script>

<template>
  <div class="flex flex-col items-center mb-8">
    <img :src="props.logoSrc" alt="Nexus" class="h-12 w-auto mb-2" />
    <p class="text-sm text-gray-700">{{ props.subtitle }}</p>
  </div>
</template>