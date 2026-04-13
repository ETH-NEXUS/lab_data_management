<script setup lang="ts">
type Props = {
  stockId: number
  productName: string
}

const props = defineProps<Props>()

const emit = defineEmits<{
  (e: 'close'): void
}>()

const { t } = useI18n()
</script>

<template>
  <header class="flex items-start justify-between gap-4 border-b border-[var(--app-border)] px-5 py-4">
    <div class="min-w-0 space-y-1">
      <p class="text-xs font-semibold tracking-[0.12em] text-slate-500 uppercase">
        {{ t('inventory.stock_drawer.title') }}
      </p>
      <h2 class="truncate text-lg font-semibold text-slate-900">
        {{ props.productName }}
      </h2>
      <p class="text-xs text-slate-600">
        {{ t('inventory.stock_drawer.subtitle', { id: props.stockId }) }}
      </p>
    </div>

    <UButton
      size="xs"
      variant="ghost"
      color="neutral"
      icon="i-heroicons-x-mark"
      :label="t('inventory.stock_drawer.actions.close')"
      @click="emit('close')"
    />
  </header>
</template>
