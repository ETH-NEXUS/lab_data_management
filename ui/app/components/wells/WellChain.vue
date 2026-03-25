<script setup lang="ts">
import { nextTick, onMounted, reactive, ref } from 'vue'
import * as vNG from 'v-network-graph'
import type { Edges, Layouts, Nodes } from 'v-network-graph'
import { useAPI } from '~/composables/useAPI'
import type { Well } from '~/types/lab'

const props = defineProps<{
  well: Well
}>()

const loading = ref<boolean>(true)
const nodes = ref<Nodes>({})
const edges = ref<Edges>({})
const layouts: Layouts = reactive({
  nodes: {},
})

const graph = ref<vNG.VNetworkGraphInstance>()
const VNetworkGraph = vNG.VNetworkGraph
const VEdgeLabel = vNG.VEdgeLabel

const nodeSize = 20

const configs = vNG.defineConfigs({
  view: {
    panEnabled: true,
    zoomEnabled: true,
  },
  node: {
    draggable: false,
    normal: {
      radius: nodeSize / 2,
      color: (node: unknown) => {
        const value = node as { root?: boolean }
        return value.root ? '#ff0000' : '#4466cc'
      },
    },
    label: { direction: 'south', color: '#333' },
  },
  edge: {
    normal: {
      color: '#aaa',
      width: 3,
    },
    margin: 4,
    marker: {
      target: {
        type: 'arrow',
        width: 4,
        height: 4,
      },
    },
  },
})

const layout = (direction: 'TB' | 'LR') => {
  if (Object.keys(nodes.value).length <= 1 || Object.keys(edges.value).length === 0) {
    return
  }

  const nodeIds = Object.keys(nodes.value)
  const out = new Map<string, string[]>()
  const indegree = new Map<string, number>()

  for (const nodeId of nodeIds) {
    out.set(nodeId, [])
    indegree.set(nodeId, 0)
  }

  for (const edge of Object.values(edges.value)) {
    out.get(edge.source)?.push(edge.target)
    indegree.set(edge.target, (indegree.get(edge.target) ?? 0) + 1)
  }

  const queue: string[] = nodeIds.filter((id) => (indegree.get(id) ?? 0) === 0)
  const level = new Map<string, number>(queue.map((id) => [id, 0]))

  while (queue.length > 0) {
    const current = queue.shift()
    if (!current) continue

    const currentLevel = level.get(current) ?? 0
    for (const target of out.get(current) ?? []) {
      const nextLevel = currentLevel + 1
      const previous = level.get(target)
      if (previous === undefined || nextLevel > previous) {
        level.set(target, nextLevel)
      }

      indegree.set(target, (indegree.get(target) ?? 1) - 1)
      if ((indegree.get(target) ?? 0) === 0) {
        queue.push(target)
      }
    }
  }

  const buckets = new Map<number, string[]>()
  for (const nodeId of nodeIds) {
    const lvl = level.get(nodeId) ?? 0
    if (!buckets.has(lvl)) buckets.set(lvl, [])
    buckets.get(lvl)?.push(nodeId)
  }

  const box: Record<string, number | undefined> = {}
  const levels = [...buckets.keys()].sort((a, b) => a - b)
  levels.forEach((lvl, levelIndex) => {
    const ids = buckets.get(lvl) ?? []
    ids.forEach((nodeId, rowIndex) => {
      const x = direction === 'LR' ? levelIndex * nodeSize * 6 : rowIndex * nodeSize * 4
      const y = direction === 'LR' ? rowIndex * nodeSize * 4 : levelIndex * nodeSize * 6
      layouts.nodes[nodeId] = { x, y }

      box.top = box.top ? Math.min(box.top, y) : y
      box.bottom = box.bottom ? Math.max(box.bottom, y) : y
      box.left = box.left ? Math.min(box.left, x) : x
      box.right = box.right ? Math.max(box.right, x) : x
    })
  })

  const graphMargin = nodeSize * 2
  const viewBox = {
    top: (box.top ?? 0) - graphMargin,
    bottom: (box.bottom ?? 0) + graphMargin,
    left: (box.left ?? 0) - graphMargin,
    right: (box.right ?? 0) + graphMargin,
  }
  graph.value?.setViewBox(viewBox)
  graph.value?.panToCenter()
  graph.value?.fitToContents()
}

onMounted(async () => {
  try {
    const { data, error } = await useAPI<{ nodes: Nodes; edges: Edges }>(`wells/${props.well.id}/chain/`, {
      method: 'GET',
    })
    if (!error.value && data.value) {
      nodes.value = data.value.nodes
      edges.value = data.value.edges
      nextTick(() => {
        layout('LR')
      })
    }
  } catch (err) {
    console.error(err)
  } finally {
    loading.value = false
  }
})
</script>

<template>
  <div class="w-full" style="height: 300px">
    <div v-if="loading" class="flex h-full items-center justify-center">
      <span class="inline-block h-8 w-8 animate-spin rounded-full border-2 border-teal-800/30 border-t-teal-800" />
    </div>
    <ClientOnly>
      <VNetworkGraph
        v-if="!loading"
        ref="graph"
        class="graph h-full w-full"
        :nodes="nodes"
        :edges="edges"
        :layouts="layouts"
        :configs="configs"
      >
        <template #edge-label="{ edge, ...slotProps }">
          <VEdgeLabel :text="edge.label" align="center" vertical-align="above" v-bind="slotProps" />
        </template>
      </VNetworkGraph>
    </ClientOnly>
  </div>
</template>
