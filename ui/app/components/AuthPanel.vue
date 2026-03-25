<script setup lang="ts">
const authStore = useAuthStore()

const username = ref('')
const password = ref('')

const authErrorStatus = useState<number | null>('authErrorStatus', () => null)

const submitDisabled = computed(() => authStore.isLoading || !username.value || !password.value)

const onLogin = async () => {
  authStore.clearError()
  const ok = await authStore.login(username.value, password.value)
  if (ok) {
    password.value = ''
  }
}
</script>

<template>

  <UCard class="w-full">
    <template #header>
      <div class="flex items-center justify-between gap-2">
        <h3 class="text-base font-semibold">Auth</h3>
        <UBadge v-if="authStore.isAuthenticated" color="success" variant="soft">signed in</UBadge>
        <UBadge v-else color="info" variant="soft">signed out</UBadge>
      </div>
    </template>

    <div class="flex flex-col gap-3">
      <div v-if="authStore.isAuthenticated" class="flex flex-col gap-2">
        <div class="text-sm">
          <div v-if="authStore.user?.username">
            <span class="font-semibold">username:</span> {{ authStore.user.username }}
          </div>
          <div v-if="authStore.user?.email"><span class="font-semibold">email:</span> {{ authStore.user.email }}</div>
        </div>

        <div class="flex flex-wrap items-center gap-2">
          <UButton size="xs" variant="soft" :loading="authStore.isLoading" @click="authStore.loadUser()">
            Refresh user
          </UButton>
          <UButton size="xs" color="red" variant="soft" :loading="authStore.isLoading" @click="authStore.logout()">
            Logout
          </UButton>
        </div>
      </div>

      <form v-else class="flex flex-col gap-2" @submit.prevent="onLogin">
        <UInput v-model="username" placeholder="Username" autocomplete="username" />
        <UInput v-model="password" type="password" placeholder="Password" autocomplete="current-password" />

        <div class="flex flex-wrap items-center gap-2">
          <UButton type="submit" size="xs" variant="soft" :loading="authStore.isLoading" :disabled="submitDisabled">
            Login
          </UButton>
          <UButton size="xs" variant="ghost" :loading="authStore.isLoading" @click="authStore.loadUser()">
            Check session
          </UButton>
        </div>
      </form>

      <UBadge v-if="authStore.error" color="red" variant="soft">{{ authStore.error }}</UBadge>
      <UBadge v-else-if="authErrorStatus" color="amber" variant="soft">Auth error: {{ authErrorStatus }}</UBadge>
    </div>
  </UCard>
</template>
