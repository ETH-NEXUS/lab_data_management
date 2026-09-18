import { fileURLToPath } from 'node:url'
import { defineConfig } from 'vitest/config'

// The unit tests run without Nuxt: they import a store or a helper directly,
// so "~" has to point at this folder the way Nuxt does it.
export default defineConfig({
  resolve: {
    alias: {
      '~': fileURLToPath(new URL('.', import.meta.url)),
    },
  },
  test: {
    environment: 'happy-dom',
    include: ['tests/**/*.test.ts'],
  },
})
