import VueExcelEditor from 'vue3-excel-editor'

/**
 * Registers the Excel editor plugin on client side.
 *
 * Registered global components:
 * - `<vue-excel-editor />`
 * - `<vue-excel-column />`
 */
export default defineNuxtPlugin((nuxtApp) => {
  nuxtApp.vueApp.use(VueExcelEditor)
})
