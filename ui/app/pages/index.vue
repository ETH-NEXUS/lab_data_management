<script setup lang="ts">
import ldmLogo from '~/assets/ldm_logo.png'
</script>

<template>
  <section class="px-4 py-10 sm:px-6 lg:px-8">
    <div class="mx-auto flex min-h-[60vh] max-w-4xl items-center justify-center">
      <div class="w-full max-w-2xl rounded-3xl border border-slate-200/80 bg-white/90 px-8 py-14 text-center shadow-sm">
        <div class="mb-8 inline-flex rounded-2xl bg-slate-50 p-4">
          <img :src="ldmLogo" alt="LDM Logo" class="h-16 w-auto object-contain sm:h-20" />
        </div>

        <p class="text-xs font-semibold tracking-[0.16em] text-blue-600 uppercase">Dashboard</p>
        <h1 class="mt-3 text-3xl font-semibold tracking-tight text-slate-900 sm:text-4xl">Lab Management</h1>
        <p class="mt-4 text-base text-slate-600 sm:text-lg">
          Open any section from the left navigation to continue your workflow.
        </p>
      </div>
    </div>
  </section>
</template>
