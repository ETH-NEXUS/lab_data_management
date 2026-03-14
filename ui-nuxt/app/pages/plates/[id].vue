<script setup lang="ts">
import { computed, ref } from 'vue'
import PlatePaneHeader from '~/components/plates/PlatePaneHeader.vue'
import ResizableSplitPane from '~/components/common/ResizableSplitPane.vue'
import { PLATE_PAGE_DEFAULT_LEFT_PERCENT, PLATE_PAGE_MAX_LEFT_PERCENT, PLATE_PAGE_MIN_LEFT_PERCENT } from '~/types/plates'
import { getPlateRouteIdLabel } from '~/utils/plates'

const route = useRoute()
const { t } = useI18n()
const leftPanelPercent = ref(PLATE_PAGE_DEFAULT_LEFT_PERCENT)
const plateIdLabel = computed(() => getPlateRouteIdLabel(route.params.id))
</script>

<template>
  <section class="h-[calc(100dvh-0.75rem)] px-3 pb-3 pt-2">
    <div class="h-full border border-black/10 bg-white/25 backdrop-blur-sm">
      <ResizableSplitPane
        v-model="leftPanelPercent"
        :min-left-percent="PLATE_PAGE_MIN_LEFT_PERCENT"
        :max-left-percent="PLATE_PAGE_MAX_LEFT_PERCENT"
        :divider-width="3"
      >
        <template #left>
          <div class="h-full p-4">
            <div class="h-full">
              <PlatePaneHeader
                :title="t('plates.page.title', { id: plateIdLabel })"
                :description="t('plates.page.left_placeholder')"
                title-class="mb-3 font-mono text-3xl font-medium text-slate-800"
              />
            </div>
          </div>
        </template>

        <template #right>
          <div class="h-full p-4">
            <div class="h-full">
              <PlatePaneHeader
                :title="t('plates.page.right_title')"
                :description="t('plates.page.right_placeholder')"
              />
            </div>
          </div>
        </template>
      </ResizableSplitPane>
    </div>
  </section>
</template>
