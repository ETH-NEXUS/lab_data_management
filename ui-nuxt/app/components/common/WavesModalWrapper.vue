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
  footerClass: 'flex w-full flex-wrap justify-start gap-2',
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
  overlay: 'bg-slate-900/55 backdrop-blur-[2px]',
  content:
    'overflow-hidden rounded-3xl border border-blue-200/80 bg-white p-0 shadow-[0_30px_90px_rgba(30,64,175,0.28)]',
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
      <section class="relative overflow-hidden">
        <div class="relative overflow-hidden bg-blue-500 px-4 pt-12 pb-8 sm:px-8 sm:pt-14 sm:pb-10">
          <div class="pointer-events-none absolute -top-8 -right-16 h-40 w-40 rounded-full border-2 border-white/25" />
          <div
            class="pointer-events-none absolute -bottom-10 -left-16 h-44 w-44 rounded-full border-2 border-orange-300/45"
          />
          <div class="pointer-events-none absolute right-8 bottom-5 h-14 w-14 rounded-full border border-white/20" />

          <header class="mx-auto max-w-2xl text-center">
            <h2 class="text-2xl leading-tight font-bold text-white sm:text-4xl">
              {{ props.title }}
            </h2>
            <p v-if="props.description" class="mt-3 text-sm text-blue-100 sm:text-base">
              {{ props.description }}
            </p>
          </header>

          <div class="relative mx-auto mt-8 max-w-2xl rounded-xl border border-slate-200 bg-white p-1 shadow-lg">
            <div :class="['mx-auto', props.bodyContainerClass]">
              <slot name="body" />
            </div>

            <div v-if="$slots.footer" :class="['mt-4 border-t border-slate-200 bg-white px-6 py-5', props.footerClass]">
              <slot name="footer" />
            </div>
          </div>
        </div>
      </section>
    </template>
  </UModal>
</template>
