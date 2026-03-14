import tailwindcss from '@tailwindcss/vite'
import { fileURLToPath } from 'node:url'

const srcDir = fileURLToPath(new URL('.', import.meta.url))
const apiBaseURL = process.env.NUXT_PUBLIC_API_URL ? `${process.env.NUXT_PUBLIC_API_URL}/api` : '/api'


// https://nuxt.com/docs/api/configuration/nuxt-config
export default defineNuxtConfig({
  srcDir: '.',
  ssr: false,
  css: ['./assets/css/main.css'],
  nitro: {
    static: true,
    devProxy: {
      '/api': `http://api:${process.env.DJANGO_PORT || '5000'}/api`,
      '/admin': `http://api:${process.env.DJANGO_PORT || '5000'}/admin`,
      '/media': `http://api:${process.env.DJANGO_PORT || '5000'}/media`,
      '/static': `http://api:${process.env.DJANGO_PORT || '5000'}/static`,
      '/docs': `http://api:${process.env.DJANGO_PORT || '5000'}/docs`,
      '/notebook': {
        target: `http://api:${process.env.DJANGO_PORT || '5000'}/notebook`,
        ws: true,
      },
    },
  },
  runtimeConfig: {
    public: {
      baseURL: process.env.NUXT_PUBLIC_API_URL || '',
    },
  },
  vite: {
    plugins: [tailwindcss()],
  },
  devtools: { enabled: true },
  modules: ['@nuxt/ui', '@nuxt/eslint', '@pinia/nuxt', '@vueuse/nuxt', '@nuxtjs/i18n'],
  i18n: {
    restructureDir: '',
    strategy: 'prefix_except_default',
    defaultLocale: 'en',
    detectBrowserLanguage: {
      useCookie: true,
      cookieKey: 'i18n_redirected',
      redirectOn: 'root',
    },
    locales: [
      {
        code: 'en',
        name: 'English',
        file: 'en.json',
      },
    ],
    langDir: 'locales',
  },
  devServer: {
    host: '0.0.0.0',
    port: parseInt(process.env.UI_PORT || '8077', 10),
  },
  alias: {
    components: `${srcDir}/components`,
    stores: `${srcDir}/stores`,
    types: `${srcDir}/types`,
    utils: `${srcDir}/utils`,
    errors: `${srcDir}/utils/errors`,
  },
})
