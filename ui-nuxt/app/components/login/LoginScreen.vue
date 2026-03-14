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
  <BackgroundPageWaves background-src="/assets/double-waves.png">
    <div class="w-full max-w-sm">
      <form @submit.prevent="onLogin">
        <h3 class="text-4xl text-center font-medium mb-10">
          {{ t('auth.login.title') }}
        </h3>

        <AuthBrandHeader
          logo-src="/assets/nexus_logo.png"
          :subtitle="t('auth.login.subtitle')"
        />

        <div
          v-if="showLoginError"
          class="mb-6 px-4 py-3 rounded-2xl bg-red-200/70 text-red-900 shadow"
        >
          <div class="font-medium">{{ t('auth.login.error_title') }}</div>
          <div class="text-sm">{{ t('auth.login.error_description') }}</div>
        </div>

        <BaseField
          v-model="username"
          :label="t('auth.login.username_label')"
          :placeholder="t('auth.login.username_placeholder')"
          autocomplete="username"
          class="mb-6"
          field-class=""
        />

        <BaseField
          v-model="password"
          :label="t('auth.login.password_label')"
          :placeholder="t('auth.login.password_placeholder')"
          :type="passwordVisible ? 'text' : 'password'"
          autocomplete="current-password"
          :has-right-slot="true"
          class="mb-6"
        >
          <template #right>
            <button
              type="button"
              class="p-1"
              @click="passwordVisible = !passwordVisible"
              :aria-label="passwordVisible ? 'Hide password' : 'Show password'"
            >
              <EyeIcon />
            </button>
          </template>
        </BaseField>

        <div class="text-right mb-10">
          <a href="#" class="inline-block text-sm underline font-medium">
            {{ t('auth.login.forgot_password') }}
          </a>
        </div>

        <BaseButton
          :label="t('auth.login.submit')"
          :on-click="onLogin"
          :loading="authStore.isLoading"
          :disabled="submitDisabled"
        />

        <div v-if="authStore.error" class="mt-4 text-center text-sm text-red-800">
          {{ authStore.error }}
        </div>
      </form>
    </div>
  </BackgroundPageWaves>
</template>
