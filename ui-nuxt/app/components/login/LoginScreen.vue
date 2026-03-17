<script setup lang="ts">
import { ref, computed } from 'vue'
import { useRoute, useRouter } from 'vue-router'
import BackgroundPageWaves from '~/components/common/BackgroundPageWaves.vue'
import AuthBrandHeader from '~/components/login/AuthHeader.vue'
import BaseField from '~/components/common/BaseField.vue'
import BaseButton from '~/components/common/BaseButton.vue'
import EyeIcon from '~/components/icons/EyeIcon.vue'

const { t } = useI18n()

const authStore = useAuthStore()
const router = useRouter()
const route = useRoute()

const username = ref('')
const password = ref('')
const showLoginError = ref(false)
const passwordVisible = ref(false)

const submitDisabled = computed(() => {
  return authStore.isLoading || !username.value.trim() || !password.value
})

const nextUrl = computed(() => {
  const next = route.query.next
  return typeof next === 'string' && next.length > 0 ? next : '/'
})

const onLogin = async () => {
  showLoginError.value = false
  authStore.clearError()

  const ok = await authStore.login(username.value.trim(), password.value)

  if (ok) {
    password.value = ''
    await router.push({ path: nextUrl.value })
  } else {
    showLoginError.value = true
  }
}
</script>

<template>
  <BackgroundPageWaves container-class="mx-auto w-full max-w-6xl px-4 py-10 sm:px-6 lg:px-8">
    <section class="relative w-full">
      <div class="pointer-events-none absolute inset-x-2 top-8 bottom-8 rounded-3xl bg-blue-200/80 sm:inset-x-6" />

      <div
        class="relative grid gap-8 rounded-3xl border border-blue-300/70 bg-blue-500 p-6 text-white shadow-[0_24px_60px_rgba(30,64,175,0.25)] sm:p-8 lg:grid-cols-[1.05fr_0.95fr] lg:p-10"
      >
        <div class="flex flex-col justify-center">
          <span
            class="mb-4 inline-flex w-fit rounded-full border border-white/30 bg-white/10 px-4 py-1 text-xs font-semibold tracking-[0.16em] uppercase"
          >
            Nexus Lab
          </span>

          <h1 class="text-3xl font-semibold tracking-tight sm:text-4xl">
            {{ t('auth.login.title') }}
          </h1>

          <p class="mt-4 max-w-lg text-base leading-relaxed text-blue-100 sm:text-lg">
            {{ t('auth.login.subtitle') }}
          </p>

          <div class="mt-8 grid gap-3 text-sm text-blue-100 sm:text-base">
            <p class="flex items-center gap-2">
              <span class="h-2 w-2 rounded-full bg-blue-100" />
              Manage projects, plates, and experiment results from one place.
            </p>
            <p class="flex items-center gap-2">
              <span class="h-2 w-2 rounded-full bg-blue-100" />
              Keep your workflow fast with a focused and clean workspace.
            </p>
          </div>
        </div>

        <div class="rounded-2xl border border-white/70 bg-white p-6 shadow-lg sm:p-8">
          <form @submit.prevent="onLogin">
            <AuthBrandHeader logo-src="/assets/nexus_logo.png" :subtitle="t('auth.login.subtitle')" />

            <div v-if="showLoginError" class="mb-5 rounded-xl border border-red-300 bg-red-50 px-4 py-3 text-red-900">
              <div class="font-medium">{{ t('auth.login.error_title') }}</div>
              <div class="text-sm">{{ t('auth.login.error_description') }}</div>
            </div>

            <BaseField
              v-model="username"
              :label="t('auth.login.username_label')"
              :placeholder="t('auth.login.username_placeholder')"
              autocomplete="username"
              class="mb-5"
            />

            <BaseField
              v-model="password"
              :label="t('auth.login.password_label')"
              :placeholder="t('auth.login.password_placeholder')"
              :type="passwordVisible ? 'text' : 'password'"
              autocomplete="current-password"
              :has-right-slot="true"
              class="mb-5"
            >
              <template #right>
                <button
                  type="button"
                  class="rounded-md p-1 text-slate-500 transition hover:bg-slate-100 hover:text-slate-800"
                  :aria-label="passwordVisible ? 'Hide password' : 'Show password'"
                  @click="passwordVisible = !passwordVisible"
                >
                  <EyeIcon />
                </button>
              </template>
            </BaseField>

            <div class="mb-7 text-right">
              <a
                href="#"
                class="inline-block text-sm font-medium text-[var(--app-accent)] hover:text-[var(--app-accent-hover)]"
              >
                {{ t('auth.login.forgot_password') }}
              </a>
            </div>

            <BaseButton
              :label="t('auth.login.submit')"
              :on-click="onLogin"
              :loading="authStore.isLoading"
              :disabled="submitDisabled"
              class-name="h-12 rounded-xl text-sm font-semibold uppercase tracking-[0.12em]"
            />

            <div v-if="authStore.error" class="mt-4 text-center text-sm text-red-800">
              {{ authStore.error }}
            </div>
          </form>
        </div>
      </div>
    </section>
  </BackgroundPageWaves>
</template>
