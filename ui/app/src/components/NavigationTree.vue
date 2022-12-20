<script setup lang="ts">
import {ref, onMounted, computed} from 'vue'
import {api} from '../boot/axios'
import {handleError} from '../helpers/errorHandling'
import {CompoundLibrary} from './models'
import {QTreeNode} from 'quasar'
import {useI18n} from 'vue-i18n'
import {useRouter} from 'vue-router'

const router = useRouter()
const {t} = useI18n()

const compoundLibraries = ref<Array<CompoundLibrary> | undefined>()

onMounted(async () => {
    try {
        const resp = await api.get('/api/compoundlibraries/')
        compoundLibraries.value = resp.data.results
    } catch (err) {
        handleError(err)
    }
})

const expanded = ref([])
const filter = ref('')
const filterRef = ref<HTMLInputElement>()

const resetFilter = () => {
    filter.value = ''
    filterRef.value?.focus()
}

const nodeHandler = (node: QTreeNode) => {
    if ('plate' in node) {
        router.push(`/plate/${node.plate.barcode}`)
    }
}

const nodes = computed<Array<QTreeNode>>(() => {
    const _compoundLibraries: QTreeNode = {
        label: t('label.compound_libraries'),
        icon: 'science',
        header: 'compounds',
        children: [],
    }
    const _projects: QTreeNode = {
        label: t('label.projects'),
        icon: 'biotech',
        header: 'projects',
        children: [],
    }
    const nodes: Array<QTreeNode> = []
    if (compoundLibraries.value) {
        for (const library of compoundLibraries.value) {
            const node: QTreeNode = {
                label: library.name,
                icon: 'science',
                children: [],
            }
            for (const plate of library.plates) {
                node.children?.push({
                    label: `${plate.barcode} (${plate.dimension})`,
                    icon: 'o_view_module',
                    handler: nodeHandler,
                    plate: plate,
                })
            }
            _compoundLibraries.children?.push(node)
        }
    }
    nodes.push(_compoundLibraries)
    nodes.push(_projects)
    return nodes
})
</script>

<template>
    <div class="q-pa-md q-gutter-sm">
        <q-input ref="filterRef" v-model="filter" :label="t('label.filter')">
            <template v-slot:append>
                <q-icon v-if="filter !== ''" name="clear" class="cursor-pointer" @click="resetFilter" />
            </template>
        </q-input>
        <q-tree :nodes="nodes" dense node-key="label" v-model:expanded="expanded" :filter="filter">
            <template v-slot:header-projects="prop">
                <q-icon :name="prop.node.icon || 'star'" size="28px" class="q-mr-sm" />
                <q-menu touch-position context-menu>
                    <q-list dense style="min-width: 100px">
                        <q-item clickable v-close-popup>
                            <q-item-section>{{ t('action.new_project') }}</q-item-section>
                        </q-item>
                        <q-separator />
                        <q-item clickable v-close-popup>
                            <q-item-section>{{ t('action.project_properties') }}</q-item-section>
                        </q-item>
                    </q-list>
                </q-menu>
                {{ prop.node.label }}
            </template>
        </q-tree>
    </div>
</template>

<style lang="sass">
.q-tree__node-header-content
  cursor: pointer
</style>
