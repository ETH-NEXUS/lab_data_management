<script setup lang="ts">
import type { PlatePaneHeaderProps } from '~/types/plates'

const props = withDefaults(defineProps<PlatePaneHeaderProps>(), {
  titleClass: 'mb-3 text-3xl font-semibold tracking-tight text-slate-900',
  descriptionClass: 'text-sm text-slate-700',
})
</script>

<template>
  <div>
    <h2 :class="props.titleClass">
      {{ props.title }}
    </h2>
    <p :class="props.descriptionClass">
      {{ props.description }}
    </p>
  </div>
</template>

