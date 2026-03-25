<script setup lang="ts">
import ldmLogo from '~/assets/ldm_logo.png'
</script>

<template>
  <section class="px-4 py-6 sm:px-6 lg:px-8">
    <div class="relative mx-auto max-w-6xl">
      <div class="pointer-events-none absolute inset-x-4 top-8 bottom-8 rounded-3xl bg-blue-200/70 sm:inset-x-8" />

      <div
        class="relative overflow-hidden rounded-3xl border border-blue-300/70 bg-blue-500 px-6 py-8 text-white shadow-[0_24px_60px_rgba(30,64,175,0.25)] sm:px-10 sm:py-10"
      >
        <div class="grid gap-8 lg:grid-cols-[auto_1fr_auto] lg:items-center">
          <div class="inline-flex w-fit rounded-2xl bg-white p-3">
            <img :src="ldmLogo" alt="LDM Logo" class="h-12 w-auto object-contain sm:h-14" />
          </div>

          <div>
            <p class="mb-2 text-xs font-semibold tracking-[0.16em] text-blue-100 uppercase">Projects</p>
            <h1 class="text-3xl font-semibold tracking-tight sm:text-4xl">Project workspace overview</h1>
            <p class="mt-4 max-w-2xl text-base text-blue-100 sm:text-lg">
              Select a project from the left tree to open its summary, edit metadata, and inspect related experiments.
            </p>
          </div>

          <div class="grid gap-3 sm:grid-cols-2 lg:grid-cols-1">
            <NuxtLink
              to="/"
              class="inline-flex items-center justify-center rounded-xl border border-white/50 bg-white/10 px-5 py-3 text-sm font-semibold tracking-[0.12em] uppercase transition hover:bg-white/20"
            >
              Dashboard
            </NuxtLink>
            <NuxtLink
              to="/experiments"
              class="inline-flex items-center justify-center rounded-xl bg-white px-5 py-3 text-sm font-semibold tracking-[0.12em] text-blue-700 uppercase transition hover:bg-blue-50"
            >
              Experiments
            </NuxtLink>
          </div>
        </div>
      </div>
    </div>
  </section>
</template>
