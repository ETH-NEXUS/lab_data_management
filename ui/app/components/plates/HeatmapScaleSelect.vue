<script setup lang="ts">
import { usePlateViewStore } from '~/stores/plateView'
import { HEATMAP_SCALES } from '~/types/plates'

// The colour scale of the heatmaps of the plate and the experiment page
const plateViewStore = usePlateViewStore()
const { t } = useI18n()
</script>

<template>
  <div>
    <label class="mb-1 block pl-1 text-sm font-medium text-slate-700">
      {{ t('plates.dynamic.controls.select_color_scale') }}
    </label>
    <select
      v-model="plateViewStore.heatmapScale"
      class="w-full cursor-pointer rounded-full border border-black/15 bg-white/70 px-4 py-2 text-sm ring-offset-0 outline-none focus:ring-2 focus:ring-lime-500"
    >
      <option v-for="scale in HEATMAP_SCALES" :key="`scale-${scale}`" :value="scale">
        {{ t(`plates.dynamic.controls.color_scales.${scale}`) }}
      </option>
    </select>
  </div>
</template>
