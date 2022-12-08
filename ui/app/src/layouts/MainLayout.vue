<script setup lang="ts">
import {ref} from 'vue'
import NavigationTree from '../components/NavigationTree.vue'
import {useRoute} from 'vue-router'

const route = useRoute()
const leftDrawerOpen = ref(false)

function toggleLeftDrawer() {
    leftDrawerOpen.value = !leftDrawerOpen.value
}
</script>

<template>
    <q-layout view="hHh Lpr lFf">
        <q-header elevated>
            <q-toolbar>
                <q-btn flat dense round icon="menu" aria-label="Menu" @click="toggleLeftDrawer" />
                <q-toolbar-title>Lab Management</q-toolbar-title>
                <div>DRAFT</div>
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
