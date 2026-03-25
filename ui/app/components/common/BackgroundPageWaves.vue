<script setup lang="ts">
type Props = {
  sectionClass?: string // Extra classes for the outer section
  containerClass?: string // Extra classes for the centered container
}

const props = withDefaults(defineProps<Props>(), {
  sectionClass: 'bg-[var(--app-bg)]',
  containerClass: 'mx-auto w-full px-6 py-10 sm:px-8',
})
</script>

<template>
  <section class="relative min-h-screen overflow-hidden" :class="props.sectionClass">
    <!-- Decorative gradient accents to keep depth without image assets. -->
    <div class="pointer-events-none absolute -top-16 -left-12 h-64 w-64 rounded-full bg-[var(--app-accent-soft)] blur-3xl" />
    <div class="pointer-events-none absolute right-0 bottom-0 h-72 w-72 rounded-full bg-white/60 blur-3xl" />
    <!-- Centered layouts -->
    <div class="relative min-h-screen flex items-center justify-center" :class="props.containerClass">
      <slot />
    </div>
  </section>
</template>
