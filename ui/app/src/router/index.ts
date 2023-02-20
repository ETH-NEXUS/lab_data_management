import {route} from 'quasar/wrappers'
import {createMemoryHistory, createRouter, createWebHashHistory, createWebHistory} from 'vue-router'

import routes from './routes'
import {useUserStore} from 'src/stores/user'

/*
 * If not building with SSR mode, you can
 * directly export the Router instantiation;
 *
 * The function below can be async too; either use
 * async/await or return a Promise which resolves
 * with the Router instance.
 */

export default route(function (/* { store, ssrContext } */) {
  const createHistory = process.env.SERVER
    ? createMemoryHistory
    : process.env.VUE_ROUTER_MODE === 'history'
    ? createWebHistory
    : createWebHashHistory

  const Router = createRouter({
    scrollBehavior: () => ({left: 0, top: 0}),
    routes,

    // Leave this as is and make changes in quasar.conf.js instead!
    // quasar.conf.js -> build -> vueRouterMode
    // quasar.conf.js -> build -> publicPath
    history: createHistory(process.env.VUE_ROUTER_BASE),
  })

  Router.beforeEach((to, from, next) => {
    const userStore = useUserStore()

    if (to.path.startsWith('/admin/')) {
      next()
    }

    if (!['/login'].includes(to.path ? to.path.toString() : '')) {
      // If we do not want to login or register we check if the JWT token exists, if not we route to the login page.
      console.debug('jwt', userStore.jwt)
      if (!userStore.jwt) {
        next({
          path: '/login',
          query: {next: encodeURI(to.fullPath)},
        })
      } else {
        next()
      }
    } else {
      if (userStore.jwt) {
        // It makes no sense to login if we are already logged in so we route back to the path we come from.
        next({
          path: from.path,
        })
      } else {
        next()
      }
    }
  })

  return Router
})
