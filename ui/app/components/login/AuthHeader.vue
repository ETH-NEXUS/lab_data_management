<script setup lang="ts">
type Props = {
  logoSrc: string
  subtitle: string
}

const props = defineProps<Props>()
</script>

<template>
  <div class="mb-6">
    <div class="inline-flex items-center gap-3 rounded-xl bg-[var(--app-surface-muted)] p-2 pr-4">
      <img :src="props.logoSrc" alt="Nexus" class="h-9 w-auto rounded-md bg-white p-1.5" />
      <span class="text-xs font-semibold tracking-[0.16em] text-[var(--app-text-muted)] uppercase">
        Lab Management
      </span>
    </div>
    <p class="mt-3 text-sm text-[var(--app-text-muted)]">{{ props.subtitle }}</p>
  </div>
</template>
