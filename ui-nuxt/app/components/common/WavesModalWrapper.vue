<script setup lang="ts">
type Props = {
  open: boolean
  title: string
  description?: string
  dismissible?: boolean
  modalClass?: string
  bodyContainerClass?: string
  footerClass?: string
}

const props = withDefaults(defineProps<Props>(), {
  description: '',
  dismissible: true,
  modalClass: 'w-full sm:max-w-3xl',
  bodyContainerClass: 'w-full px-6 py-6',
  footerClass: 'flex w-full justify-start gap-2 px-6 py-5',
})

const emit = defineEmits<{
  (e: 'update:open', value: boolean): void
}>()

/**
 * Shared modal UI shell used by project and experiment modals.
 *
 * Accepted props example:
 * - `{
 *     open: true,
 *     title: 'Create Project',
 *     description: 'Choose Harvest project or custom name',
 *     modalClass: 'w-full sm:max-w-4xl'
 *   }`
 *
 * Slot usage example:
 * - `<template #body>...</template>`
 * - `<template #footer>...</template>`
 *
 * Important:
 * - Do not set `relative`/`fixed`/`absolute` here because Nuxt UI controls
 *   modal centering with its own content positioning classes.
 */
const modalUi = {
  overlay: 'bg-slate-900/45 backdrop-blur-[2px]',
  content:
    'overflow-hidden rounded-3xl border border-blue-300/70 bg-white p-0 shadow-[0_24px_80px_rgba(30,64,175,0.25)]',
}
</script>

<template>
  <UModal
    :open="props.open"
    :dismissible="props.dismissible"
    :class="props.modalClass"
    :ui="modalUi"
    @update:open="emit('update:open', $event)"
  >
    <template #body>
      <section class="relative overflow-hidden px-4 py-4 sm:px-5 sm:py-5">
        <div class="pointer-events-none absolute inset-x-3 top-6 bottom-6 rounded-3xl bg-blue-200/65" />

        <div class="relative space-y-4">
          <header class="rounded-2xl border border-blue-300/70 bg-blue-500 px-5 py-5 text-left text-white sm:px-6">
            <p class="text-[11px] font-semibold tracking-[0.14em] text-blue-100 uppercase">Nexus Modal</p>
            <h2 class="mt-1 text-2xl leading-tight font-semibold sm:text-3xl">
              {{ props.title }}
            </h2>
            <p v-if="props.description" class="mt-2 text-sm text-blue-100 sm:text-base">
              {{ props.description }}
            </p>
          </header>

          <div class="rounded-2xl border border-[var(--app-border)] bg-white shadow-sm">
            <div :class="['mx-auto', props.bodyContainerClass]">
              <slot name="body" />
            </div>
          </div>
        </div>
      </section>
    </template>

    <template v-if="$slots.footer" #footer>
      <div :class="['border-t border-[var(--app-border)] bg-white', props.footerClass]">
        <slot name="footer" />
      </div>
    </template>
  </UModal>
</template>
