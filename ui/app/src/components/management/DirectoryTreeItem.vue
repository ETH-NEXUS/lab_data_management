<script setup lang="ts">
/*

 This component is not in used at the moment, but we might need it later

*/
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
