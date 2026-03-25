<script setup lang="ts">
import {defineProps, PropType, onMounted, ref, nextTick, reactive} from 'vue'
import {Well} from '../models'
import * as vNG from 'v-network-graph'
import {Nodes, Edges, Layouts} from 'v-network-graph'
import {api} from 'boot/axios'
import {handleError} from 'src/helpers/errorHandling'
// At the moment we need to base on this...
import dagre from 'dagre/dist/dagre.min.js'

const props = defineProps({
  well: {
    type: Object as PropType<Well>,
    required: true,
  },
})

const loading = ref<boolean>(true)
const nodes = ref<Nodes>({})
const edges = ref<Edges>({})
const layouts: Layouts = reactive({
  nodes: {},
})

const graph = ref<vNG.VNetworkGraphInstance>()

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
      color: n => (n.root ? '#ff0000' : '#4466cc'),
    },
    label: {direction: 'south', color: '#333'},
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
  if (Object.keys(nodes.value).length <= 1 || Object.keys(edges.value).length == 0) {
    return
  }

  // convert graph
  // ref: https://github.com/dagrejs/dagre/wiki
  const g = new dagre.graphlib.Graph()
  // Set an object for the graph label
  g.setGraph({
    rankdir: direction,
    nodesep: nodeSize * 2,
    edgesep: nodeSize,
    ranksep: nodeSize * 6,
  })
  // Default to assigning a new object as a label for each new edge.
  g.setDefaultEdgeLabel(() => ({}))

  // Add nodes to the graph. The first argument is the node id. The second is
  // metadata about the node. In this case we're going to add labels to each of
  // our nodes.
  Object.entries(nodes.value).forEach(([nodeId, node]) => {
    g.setNode(nodeId, {label: node.name, width: nodeSize, height: nodeSize})
  })

  // Add edges to the graph.
  Object.values(edges.value).forEach(edge => {
    g.setEdge(edge.source, edge.target)
  })

  dagre.layout(g)

  const box: Record<string, number | undefined> = {}
  g.nodes().forEach((nodeId: string) => {
    // update node position
    const x = g.node(nodeId).x
    const y = g.node(nodeId).y
    layouts.nodes[nodeId] = {x, y}

    // calculate bounding box size
    box.top = box.top ? Math.min(box.top, y) : y
    box.bottom = box.bottom ? Math.max(box.bottom, y) : y
    box.left = box.left ? Math.min(box.left, x) : x
    box.right = box.right ? Math.max(box.right, x) : x
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

// const updateLayout = (direction: 'TB' | 'LR') => {
//   // Animates the movement of an element.
//   graph.value?.transitionWhile(() => {
//     layout(direction)
//   })
// }

onMounted(async () => {
  try {
    const resp = await api.get(`/api/wells/${props.well.id}/chain/`)
    nodes.value = resp.data.nodes
    edges.value = resp.data.edges
    nextTick(() => {
      layout('LR')
    })
  } catch (err) {
    handleError(err)
  } finally {
    loading.value = false
  }
})
</script>

<template>
  <div class="full-width fit row wrap justify-center items-start content-start" style="height: 300px">
    <q-spinner-dots v-if="loading" color="primary" size="2em" />
    <v-network-graph
      v-else
      ref="graph"
      class="graph"
      :nodes="nodes"
      :edges="edges"
      :layouts="layouts"
      :configs="configs">
      <template #edge-label="{edge, ...slotProps}">
        <v-edge-label :text="edge.label" align="center" vertical-align="above" v-bind="slotProps" />
      </template>
    </v-network-graph>
  </div>
</template>
