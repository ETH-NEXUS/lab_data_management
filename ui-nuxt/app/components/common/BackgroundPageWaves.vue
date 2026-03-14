<script setup lang="ts">
type Props = {

  backgroundSrc: string // Background image URL (e.g. "assets/double-waves.png")
  sectionClass?: string // Extra classes for the outer section
  containerClass?: string // Extra classes for the centered container
}

const props = withDefaults(defineProps<Props>(), {
  sectionClass: 'bg-blue-50',
  containerClass: 'container mx-auto px-4'
})
</script>

<template>
  <section class="relative min-h-screen overflow-hidden" :class="props.sectionClass">
    <!-- Background -->
    <img
      class="absolute top-0 left-0 h-full w-full object-cover pointer-events-none select-none"
      :src="props.backgroundSrc"
      alt=""
    />
    <!-- Centered layouts -->
    <div class="relative min-h-screen flex items-center justify-center" :class="props.containerClass">
      <slot />
    </div>
  </section>
</template>