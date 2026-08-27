<script setup lang="ts">
type InventoryActionCardItem = {
  id: string
  title: string
  description: string
  icon: string
}

type Props = {
  item: InventoryActionCardItem
  size?: 'quick' | 'view'
}

const props = withDefaults(defineProps<Props>(), {
  size: 'quick',
})

const emit = defineEmits<{
  (e: 'select', actionId: string): void
}>()

/**
 * Emits the selected action id when one inventory card is clicked.
 *
 * Accepted data example:
 * - `props.item = { id: 'search', title: 'Search', description: 'Find entries', icon: 'i-heroicons-magnifying-glass' }`
 *
 * Returned data example:
 * - emits `select` with `'search'`
 */
const onSelect = (): void => {
  emit('select', props.item.id)
}
</script>

<template>
  <button
    type="button"
    class="inventory-action-card group"
    :class="[props.size === 'quick' ? 'inventory-quick-action-card' : 'inventory-view-action-card']"
    @click="onSelect"
  >
    <div class="mb-4 flex items-center justify-between gap-3">
      <span class="inventory-icon-chip">
        <UIcon :name="props.item.icon" class="size-5" />
      </span>

      <UIcon
        name="i-heroicons-arrow-right"
        class="size-4 text-slate-400 transition-transform duration-200 group-hover:translate-x-0.5 group-hover:text-blue-600"
      />
    </div>

    <p class="inventory-action-card-title">{{ props.item.title }}</p>
    <p class="inventory-action-card-description">{{ props.item.description }}</p>
  </button>
</template>
