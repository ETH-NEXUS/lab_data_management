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
  bodyContainerClass: 'w-full px-8 pt-10',
  footerClass: 'flex w-full justify-start gap-2 px-8 pb-8',
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
  content:
    'overflow-hidden rounded-2xl bg-[url(/assets/double-waves.png)] bg-cover bg-center bg-no-repeat p-0 shadow-md',
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
      <section class="relative overflow-hidden pb-8 md:pb-10">
        <div :class="props.bodyContainerClass">
          <header class="mb-8 text-left">
            <h2 class="text-3xl leading-tight font-semibold text-teal-950 sm:text-4xl">
              {{ props.title }}
            </h2>
            <p v-if="props.description" class="mt-2 text-sm text-slate-700 sm:text-base">
              {{ props.description }}
            </p>
          </header>

          <slot name="body" />
        </div>
      </section>
    </template>

    <template v-if="$slots.footer" #footer>
      <div :class="props.footerClass">
        <slot name="footer" />
      </div>
    </template>
  </UModal>
</template>
