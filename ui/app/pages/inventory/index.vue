<script setup lang="ts">
/**
 * Inventory landing page.
 *
 * For now this is only a styled entry page based on the lab team's sketch.
 * The buttons are placeholders and will later be connected to real views.
 */

const showPlaceholder = (label: string): void => {
  window.alert(`${label}: to be implemented`)
}
</script>

<template>
  <section class="px-4 py-6 sm:px-6 lg:px-8">
    <div class="relative mx-auto max-w-6xl">
      <div class="pointer-events-none absolute inset-x-4 top-8 bottom-8 rounded-3xl bg-blue-200/70 sm:inset-x-8" />

      <div
        class="relative overflow-hidden rounded-3xl border border-blue-300/70 bg-blue-500 px-6 py-8 text-white shadow-[0_24px_60px_rgba(30,64,175,0.25)] sm:px-10 sm:py-10"
      >
        <div class="space-y-8">
          <div class="flex flex-col gap-3">
            <p class="text-xs font-semibold tracking-[0.16em] text-blue-100 uppercase">Laboratory inventory</p>

            <h1 class="text-3xl font-semibold tracking-tight sm:text-4xl">Inventory</h1>

            <p class="max-w-3xl text-base text-blue-100 sm:text-lg">
              A central place for stock, devices, reserved materials, recent activity, and all the items that currently
              live in the old spreadsheet workflow.
            </p>
          </div>

          <div class="rounded-3xl border border-white/20 bg-white/95 p-5 text-slate-900 shadow-inner sm:p-6">
            <div class="flex flex-col gap-8">
              <!-- Top action row -->
              <div class="flex flex-col gap-6 xl:flex-row xl:items-start xl:justify-between">
                <div class="flex flex-wrap items-center gap-3">
                  <button class="inventory-pill-button" type="button" @click="showPlaceholder('All items')">All</button>

                  <button
                    class="inventory-pill-button inventory-pill-button--icon"
                    type="button"
                    @click="showPlaceholder('Search')"
                  >
                    <span>Search</span>
                    <span class="text-lg leading-none">⌕</span>
                  </button>

                  <button
                    class="inventory-pill-button inventory-pill-button--icon"
                    type="button"
                    @click="showPlaceholder('Advanced search')"
                  >
                    <span>Advanced</span>
                    <span class="text-lg leading-none">⌕</span>
                  </button>

                  <button class="inventory-pill-button" type="button" @click="showPlaceholder('Order')">Order</button>
                </div>

                <div class="flex flex-wrap items-start gap-5">
                  <div class="flex flex-col items-center gap-2">
                    <p class="text-sm font-medium text-slate-700">Add new item</p>

                    <button
                      class="inventory-round-button"
                      type="button"
                      aria-label="Add new item"
                      @click="showPlaceholder('Add new item')"
                    >
                      <span class="text-3xl leading-none font-semibold">+</span>
                    </button>
                  </div>

                  <div class="flex flex-col items-center gap-2">
                    <button
                      class="inventory-round-button"
                      type="button"
                      aria-label="Favorite items"
                      @click="showPlaceholder('Favorite items')"
                    >
                      <span class="text-2xl leading-none">★</span>
                    </button>

                    <p class="max-w-[90px] text-center text-sm font-medium text-slate-700">Favorit Items</p>
                  </div>
                </div>
              </div>

              <!-- Main action grid -->
              <div class="grid gap-4 sm:grid-cols-2 xl:max-w-4xl">
                <button class="inventory-card-button" type="button" @click="showPlaceholder('Expired items')">
                  Expired Items
                </button>

                <button class="inventory-card-button" type="button" @click="showPlaceholder('Devices')">Devices</button>

                <button class="inventory-card-button" type="button" @click="showPlaceholder('Recent activities')">
                  Recent activities
                </button>

                <button class="inventory-card-button" type="button" @click="showPlaceholder('Recently linked items')">
                  Recently linked Items
                </button>

                <button class="inventory-card-button" type="button" @click="showPlaceholder('Reserved items')">
                  Reserved Items
                </button>

                <button class="inventory-card-button" type="button" @click="showPlaceholder('Archived items')">
                  Archived Items
                </button>
              </div>

              <div class="rounded-2xl border border-blue-100 bg-blue-50 px-4 py-4 text-sm text-slate-700">
                This is a first visual shell for the inventory page. The next step is to connect these entry points to
                the real stock, material, order, and usage components.
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  </section>
</template>

<style scoped>
.inventory-pill-button {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  gap: 10px;
  min-height: 52px;
  padding: 0 18px;
  border: 2px solid rgb(15 23 42 / 0.9);
  border-radius: 14px;
  background: white;
  color: rgb(15 23 42);
  font-size: 1.5rem;
  font-weight: 500;
  line-height: 1;
  transition:
    transform 0.15s ease,
    box-shadow 0.15s ease,
    background-color 0.15s ease;
}

.inventory-pill-button:hover {
  background: rgb(239 246 255);
  box-shadow: 0 10px 20px rgb(30 64 175 / 0.12);
  transform: translateY(-1px);
}

.inventory-pill-button:active {
  transform: translateY(0);
}

.inventory-pill-button--icon {
  justify-content: space-between;
}

.inventory-round-button {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  width: 62px;
  height: 62px;
  border: 2px solid rgb(15 23 42 / 0.9);
  border-radius: 9999px;
  background: white;
  color: rgb(15 23 42);
  transition:
    transform 0.15s ease,
    box-shadow 0.15s ease,
    background-color 0.15s ease;
}

.inventory-round-button:hover {
  background: rgb(239 246 255);
  box-shadow: 0 10px 20px rgb(30 64 175 / 0.12);
  transform: translateY(-1px);
}

.inventory-round-button:active {
  transform: translateY(0);
}

.inventory-card-button {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  min-height: 86px;
  padding: 16px 20px;
  border: 2px solid rgb(15 23 42 / 0.9);
  border-radius: 16px;
  background: white;
  color: rgb(15 23 42);
  font-size: 1.1rem;
  font-weight: 500;
  text-align: center;
  transition:
    transform 0.15s ease,
    box-shadow 0.15s ease,
    background-color 0.15s ease;
}

.inventory-card-button:hover {
  background: rgb(239 246 255);
  box-shadow: 0 12px 24px rgb(30 64 175 / 0.12);
  transform: translateY(-1px);
}

.inventory-card-button:active {
  transform: translateY(0);
}
</style>
