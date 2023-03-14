<script setup lang="ts">
import {ref, onMounted} from 'vue'
import NavigationTree from '../components/NavigationTree.vue'
import {useRoute, useRouter} from 'vue-router'
import {useUserStore} from 'src/stores/user'
import {useQuasar} from 'quasar'
import {useI18n} from 'vue-i18n'
import {api} from 'src/boot/axios'

const $q = useQuasar()
const route = useRoute()
const router = useRouter()
const userStore = useUserStore()
const leftDrawerOpen = ref(false)

const {t} = useI18n()

const toggleLeftDrawer = () => {
  leftDrawerOpen.value = !leftDrawerOpen.value
}

const logout = () => {
  userStore.removeToken()
  $q.notify({
    type: 'positive',
    message: t('message.logged_out'),
  })
  router.push({path: '/login'})
}

const openNotebook = () => {
  window.open(`${window.location.protocol}//${window.location.host}/notebook/`, '_blank', 'noreferrer')
}

const version = ref<string>('N/A')

onMounted(async () => {
  try {
    const resp = await api.get('/api/version/')
    version.value = resp.data.version
  } catch (err) {
    console.error(err)
  }
})
</script>

<template>
  <q-layout view="hHh Lpr lFf">
    <q-header elevated>
      <q-toolbar>
        <q-btn flat dense round icon="menu" aria-label="Menu" @click="toggleLeftDrawer" />
        <q-toolbar-title>Lab Data Management</q-toolbar-title>
        <span v-if="userStore.user">
          <strong>{{ userStore.user.first_name }} {{ userStore.user.last_name }}</strong>
        </span>
        <q-btn v-if="userStore.jwt" flat round dense icon="person">
          <q-menu fit>
            <q-list style="min-width: 100px">
              <q-item v-ripple>
                <q-item-section avatar>
                  <q-icon name="o_pin" />
                </q-item-section>
                <q-item-section>{{ version }}</q-item-section>
              </q-item>
              <q-separator />
              <q-item v-ripple clickable @click="openNotebook">
                <q-item-section avatar>
                  <q-icon name="o_edit_note" />
                </q-item-section>
                <q-item-section>{{ t('label.notebook') }}</q-item-section>
              </q-item>
              <q-item v-ripple clickable @click="logout">
                <q-item-section avatar>
                  <q-icon name="logout" />
                </q-item-section>
                <q-item-section>{{ t('label.logout') }}</q-item-section>
              </q-item>
            </q-list>
          </q-menu>
        </q-btn>
      </q-toolbar>
    </q-header>

    <q-drawer v-model="leftDrawerOpen" show-if-above bordered behavior="desktop" :width="400">
      <navigation-tree />
    </q-drawer>

    <q-page-container>
      <router-view :key="route.fullPath" />
    </q-page-container>
  </q-layout>
</template>
