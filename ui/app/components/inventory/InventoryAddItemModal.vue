<script setup lang="ts">
type Props = {
  open: boolean
}

const props = defineProps<Props>()

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
}>()

const closeModal = (): void => {
  emit('update:open', false)
}
</script>

<template>
  <UModal
    :open="props.open"
    title="Add new item"
    description="Choose which add-item flow you want to run."
    class="w-full sm:max-w-2xl"
    :ui="{ content: 'rounded-2xl bg-white shadow-md' }"
    @update:open="(isOpen) => emit('update:open', isOpen)"
  >
    <template #body>
      <div class="space-y-4 p-6">
        <p class="text-sm text-slate-600">Options</p>

        <div class="space-y-3">
          <div class="rounded-lg border border-slate-200 bg-slate-50 px-4 py-3">
            <p class="text-sm font-semibold text-slate-800">Add stock for existing material</p>
            <p class="mt-1 text-sm text-slate-600">Create one new row in the items table from existing catalog data.</p>
          </div>

          <div class="rounded-lg border border-slate-200 bg-slate-50 px-4 py-3">
            <p class="text-sm font-semibold text-slate-800">Create new material and stock</p>
            <p class="mt-1 text-sm text-slate-600">
              Combined creation flow will be wired in the next implementation step.
            </p>
          </div>
        </div>
      </div>
    </template>

    <template #footer>
      <div class="flex w-full justify-end gap-2 px-6 pb-6">
        <UButton variant="ghost" color="neutral" label="Cancel" @click="closeModal" />
        <UButton color="primary" label="Continue" disabled />
      </div>
    </template>
  </UModal>
</template>
