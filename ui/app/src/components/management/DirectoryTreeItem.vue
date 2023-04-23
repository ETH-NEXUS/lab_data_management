<script setup lang="ts">
import {FileSystemItem} from '../models'

defineProps<{
  item: FileSystemItem
}>()
</script>

<template>
  <li>
    {{ item.name }}
    <div v-if="item.type === 'directory'">
      <ul v-for="child in item.children" :key="child.name">
        <directory-tree-item :item="child"></directory-tree-item>
      </ul>
    </div>
  </li>
</template>
