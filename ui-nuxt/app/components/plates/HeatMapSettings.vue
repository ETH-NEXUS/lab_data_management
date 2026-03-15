<script setup lang="ts">
import { usePlateViewStore } from '~/stores/plateView'

const { t } = useI18n()
const plateViewStore = usePlateViewStore()
</script>

<template>
  <div class="mb-4 flex flex-wrap items-center gap-6">
    <label class="inline-flex cursor-pointer items-center gap-2 text-sm text-slate-700">
      <input v-model="plateViewStore.smallerMapView" type="checkbox" class="cursor-pointer" />
      <span>{{ t('plates.dynamic.settings.smaller_map_view') }}</span>
    </label>

    <label class="inline-flex cursor-pointer items-center gap-2 text-sm text-slate-700">
      <input v-model="plateViewStore.squareCompoundType" type="checkbox" class="cursor-pointer" />
      <span>{{ t('plates.dynamic.settings.show_wells_as_squares') }}</span>
    </label>
  </div>
</template>
