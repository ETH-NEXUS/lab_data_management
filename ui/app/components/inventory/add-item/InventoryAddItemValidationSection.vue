<script setup lang="ts">
type Props = {
  validationMessages: string[]
}

const props = defineProps<Props>()
</script>

<template>
  <div class="rounded-lg border border-slate-200 bg-slate-50 px-3 py-2 sm:col-span-2">
    <p class="text-xs font-semibold tracking-[0.06em] text-slate-600 uppercase">Validation</p>
    <ul v-if="props.validationMessages.length > 0" class="mt-1 list-disc space-y-0.5 pl-5 text-xs text-amber-700">
      <li v-for="validationMessage in props.validationMessages" :key="validationMessage">
        {{ validationMessage }}
      </li>
    </ul>
    <p v-else class="mt-1 text-xs text-emerald-700">Form data is valid.</p>
  </div>
</template>
