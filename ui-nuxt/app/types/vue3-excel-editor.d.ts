import type { DefineComponent } from 'vue'

declare module 'vue' {
  export interface GlobalComponents {
    VueExcelEditor: DefineComponent<Record<string, unknown>, Record<string, unknown>, unknown>
    VueExcelColumn: DefineComponent<Record<string, unknown>, Record<string, unknown>, unknown>
    'vue-excel-editor': DefineComponent<Record<string, unknown>, Record<string, unknown>, unknown>
    'vue-excel-column': DefineComponent<Record<string, unknown>, Record<string, unknown>, unknown>
  }
}

export {}
