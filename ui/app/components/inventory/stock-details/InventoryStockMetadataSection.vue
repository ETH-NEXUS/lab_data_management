<script setup lang="ts">
type MetadataField = {
  label: string
  value: string
}

type Props = {
  isExpanded: boolean
  fields: MetadataField[]
}

const props = defineProps<Props>()

const emit = defineEmits<{
  (e: 'toggle'): void
}>()

const { t } = useI18n()
</script>

<template>
  <section class="space-y-3 rounded-xl border border-slate-200 bg-white p-4">
    <div class="flex items-center justify-between gap-3">
      <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
        {{ t('inventory.stock_drawer.sections.metadata') }}
      </p>

      <UButton
        size="xs"
        variant="ghost"
        color="neutral"
        :label="
          props.isExpanded
            ? t('inventory.stock_drawer.metadata_toggle.hide')
            : t('inventory.stock_drawer.metadata_toggle.show')
        "
        @click="emit('toggle')"
      />
    </div>

    <div v-if="props.isExpanded" class="grid gap-2 sm:grid-cols-2">
      <div
        v-for="field in props.fields"
        :key="field.label"
        class="rounded-lg border border-slate-200 bg-slate-50 px-3 py-2"
      >
        <p class="text-[11px] font-semibold tracking-[0.06em] text-slate-500 uppercase">{{ field.label }}</p>
        <p class="mt-1 text-sm text-slate-800">{{ field.value }}</p>
      </div>
    </div>
  </section>
</template>
