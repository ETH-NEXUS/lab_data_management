<script setup lang="ts">
type MaterialIdentityField = {
  label: string
  value: string
  wide?: boolean
}

type Props = {
  isMaterialLoading: boolean
  fields: MaterialIdentityField[]
}

const props = defineProps<Props>()

const { t } = useI18n()
</script>

<template>
  <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
    <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
      {{ t('inventory.stock_drawer.sections.material_identity') }}
    </p>

    <p v-if="props.isMaterialLoading" class="text-sm text-slate-500">
      {{ t('inventory.stock_workspace.loading') }}
    </p>

    <div v-else class="grid gap-2 sm:grid-cols-2">
      <div
        v-for="field in props.fields"
        :key="field.label"
        :class="['rounded-lg border border-slate-200 bg-slate-50 px-3 py-2', field.wide ? 'sm:col-span-2' : '']"
      >
        <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">{{ field.label }}</p>
        <p class="mt-1 text-sm text-slate-800">{{ field.value }}</p>
      </div>
    </div>
  </section>
</template>
