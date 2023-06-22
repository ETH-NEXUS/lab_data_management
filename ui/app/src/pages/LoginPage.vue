<script setup lang="ts">
import {ref, onMounted} from 'vue'
import {useRoute, useRouter} from 'vue-router'
import {useUserStore} from 'src/stores/user'
import {useI18n} from 'vue-i18n'
import {useQuasar} from 'quasar'

const userStore = useUserStore()
const router = useRouter()
const route = useRoute()
const {t} = useI18n()
const $q = useQuasar()

const username = ref<string | null>(null)
const password = ref<string | null>(null)

const error = ref(false)

onMounted(() => {
  // api.get('/api-auth/cookie/')
})

const login = async () => {
  error.value = false
  try {
    if (username.value && password.value) {
      // await userStore.obtainToken({
      //   username: username.value,
      //   password: password.value,
      // })
      await userStore.sessionLogin({
        username: username.value,
        password: password.value,
      })
    }
  } catch (err) {
    console.error(err)
  }

  // if(userStore.jwt) {
  if (userStore.authenticated) {
    $q.notify({
      type: 'positive',
      message: t('message.successfully_logged_in'),
    })
    await router.push({path: (route.query.next as string) || '/'})
  } else {
    error.value = true
  }
}
</script>

<template>
  <div class="row justify-center">
    <div class="q-pa-md col-lg-4 col-md-6 col-sm-12 col-xs-12">
      <div class="row justify-center q-mt-lg">
        <img src="../assets/nexus_logo.png" />
        <br />
        <h3>Lab Data Management</h3>
      </div>
      <q-card v-if="error" class="bg-red-5 text-black text-bold q-mb-md">
        <q-card-section>
          {{ t('error.cannot_login') }}
        </q-card-section>
      </q-card>
      <q-card class="q-pa-sm">
        <q-form class="q-gutter-md" @submit="login">
          <q-input
            v-model="username"
            filled
            type="text"
            :label="t('label.username')"
            :hint="t('hint.username')"
            autocomplete="current-username"
            lazy-rules
            :rules="[val => (val && val.length > 0) || t('value_error.email')]" />
          <q-input
            v-model="password"
            filled
            type="password"
            :label="t('label.password')"
            :hint="t('hint.password')"
            autocomplete="current-password"
            lazy-rules
            :rules="[val => (val && val.length > 0) || t('value_error.password')]" />
          <div class="row justify-end">
            <div class="col-12">
              <q-btn class="float-right" :label="t('label.login')" type="submit" color="primary" />
            </div>
          </div>
        </q-form>
      </q-card>
    </div>
  </div>
</template>

<script lang="ts">
export default {
  name: 'LoginPage',
}
</script>

<style scoped></style>
