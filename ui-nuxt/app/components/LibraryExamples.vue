<script setup lang="ts">
import { useQuery } from '@tanstack/vue-query'

const authStore = useAuthStore()

const { t } = useI18n()

const route = useRoute()
const router = useRouter()

interface Book {
  id: number | string
  title: string
}

interface BooksResponse {
  results?: Book[]
}

const booksQuery = useQuery({
  queryKey: ['books'],
  queryFn: async () => {
    const response = await $fetch<BooksResponse>('/api/v1/books/', {
      method: 'GET',
      credentials: 'include',
    })
    return response?.results ?? []
  },
  enabled: false,
})

const bookCount = computed(() => booksQuery.data.value?.length ?? 0)
</script>

<template>
  <UCard class="w-full">
    <template #header>
      <div class="flex items-center justify-between gap-2">
        <h3 class="text-base font-semibold">{{ t('examples.library.title') }}</h3>
        <UBadge color="primary" variant="soft">vue-query</UBadge>
      </div>
    </template>

    <div class="flex flex-col gap-3">
      <div class="flex flex-wrap items-center gap-2">
        <UButton
          size="xs"
          variant="soft"
          icon="i-heroicons-arrow-path"
          :loading="booksQuery.isFetching.value"
          :disabled="!authStore.isAuthenticated"
          @click="booksQuery.refetch()"
        >
          {{ t('examples.library.refetch_books') }}
        </UButton>
        <UBadge v-if="!authStore.isAuthenticated" color="amber" variant="soft">
          {{ t('examples.library.status_login_required') }}
        </UBadge>
        <UBadge v-else-if="booksQuery.isError.value" color="red" variant="soft">{{
          t('examples.library.status_error')
        }}</UBadge>
        <UBadge v-else-if="booksQuery.isFetching.value" color="amber" variant="soft">{{
          t('examples.library.status_loading')
        }}</UBadge>
        <UBadge v-else color="green" variant="soft">{{
          t('examples.library.books_count', { count: bookCount })
        }}</UBadge>
      </div>

      <div class="flex flex-wrap items-center gap-2">
        <UBadge color="neutral" variant="soft">{{
          t('examples.library.route_label', { route: route.fullPath })
        }}</UBadge>
        <UButton size="xs" variant="ghost" @click="router.replace('/')">{{ t('examples.library.go_home') }}</UButton>
      </div>

      <ul v-if="booksQuery.data.value?.length" class="list-disc pl-5">
        <li v-for="book in booksQuery.data.value.slice(0, 3)" :key="book.id">
          {{ book.title }}
        </li>
      </ul>
    </div>
  </UCard>
</template>
