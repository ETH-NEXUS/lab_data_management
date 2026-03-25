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
      {
        path: '/management',
        component: () => import('pages/ManagementPage.vue'),
      },
      {
        path: '/messages',
        component: () => import('pages/MessagesPage.vue'),
      },
      {
        path: '/add_data/:experiment_id',
        component: () => import('pages/AddDataPage.vue'),
      },

      // await router.push(`/control-plates/${project.id}`)

      {
        path: '/control-plates/:project_id',
        component: () => import('pages/SelectControlLayoutPage.vue'),
      },
    ],
  },
  {
    path: '/notebook/:catchAll(.*)*',
    component: () => import('pages/ReloadPage.vue'),
  },
  // Always leave this as last one,
  // but you can also remove it
  {
    path: '/:catchAll(.*)*',
    component: () => import('pages/ErrorNotFound.vue'),
  },
]

export default routes
