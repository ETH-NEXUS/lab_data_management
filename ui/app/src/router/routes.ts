import {RouteRecordRaw} from 'vue-router'

const routes: RouteRecordRaw[] = [
  {
    path: '/login',
    component: () => import('pages/LoginPage.vue'),
  },
  {
    path: '/',
    component: () => import('layouts/MainLayout.vue'),
    children: [
      {
        path: '',
        component: () => import('pages/IndexPage.vue'),
      },
      {
        path: 'plate/:barcode',
        component: () => import('pages/PlatePage.vue'),
      },
      {
        path: 'project/:project',
        component: () => import('pages/ProjectPage.vue'),
      },
      {
        path: 'project/:project/experiment/:experiment',
        component: () => import('pages/ExperimentPage.vue'),
      },
    ],
  },

  // Always leave this as last one,
  // but you can also remove it
  {
    path: '/:catchAll(.*)*',
    component: () => import('pages/ErrorNotFound.vue'),
  },
]

export default routes
