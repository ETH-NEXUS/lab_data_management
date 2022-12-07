import {boot} from 'quasar/wrappers'
import axios, {AxiosInstance} from 'axios'
import {Notify, LoadingBar} from 'quasar'

declare module '@vue/runtime-core' {
    interface ComponentCustomProperties {
        $axios: AxiosInstance
    }
}

// Be careful when using SSR for cross-request state pollution
// due to creating a Singleton instance here;
// If any client changes this (global) instance, it might be a
// good idea to move this instance creation inside of the
// "export default () => {}" function below (which runs individually
// for each client)
const api = axios.create({baseURL: import.meta.env.VITE_APP_BACKEND_URL})

api.interceptors.request.use(
    config => {
        // const userStore = useUserStore()
        // const token = userStore.jwt
        // if (token) {
        //     if (config.headers) {
        //         config.headers.Authorization = `JWT ${token}`
        //     }
        // }
        // Add a trailing slash to the url
        // if (!config.url?.endsWith('/')) {
        //     config.url += '/'
        // }
        LoadingBar.start()
        return config
    },
    error => {
        console.error('AXIOS request error:', error.response)
        Notify.create({
            message: 'Error',
            caption: 'Request error',
            icon: 'warning',
            color: 'warning',
        })
        return Promise.reject(error)
    }
)

api.interceptors.response.use(
    response => {
        LoadingBar.stop()
        return response
    },
    async error => {
        // const userStore = useUserStore()
        LoadingBar.stop()
        // const originalConfig = error.config;
        // if (
        //     error.response.status === 401 &&
        //     !originalConfig._retry &&
        //     originalConfig.url !== '/auth/refresh/' &&
        //     userStore.refreshJwt
        // ) {
        //     // In case the token has expired we try to refresh the token
        //     originalConfig._retry = true
        //     await userStore.doRefreshJwt()
        //     // Resend the request with the refreshed token
        //     return api(originalConfig)
        // } else if (error.response.status === 401 && originalConfig._retry) {
        //     // If the token cannot be refreshed we logout and route to the login page
        //     await userStore.removeToken()
        //     const router = useRouter()
        //     await router.push({path: '/login'})
        // } else if (error.response) {
        if (error.response) {
            if (!error.response.data.hidden) {
                Notify.create({
                    message: `${error.response.data.detail || error.response.data.non_field_errors}`,
                    multiLine: true,
                    caption: `${error.response.config.url}: ${error.response.status}: ${error.response.statusText}`,
                    icon: 'warning',
                    color: 'negative',
                })
            }
        } else {
            Notify.create({
                message: 'Response error',
                multiLine: true,
                caption: `${error}`,
                icon: 'warning',
                color: 'negative',
            })
        }
        console.error('AXIOS response error:', error)
        return Promise.reject(error)
    }
)

export default boot(({app}) => {
    // for use inside Vue files (Options API) through this.$axios and this.$api

    app.config.globalProperties.$axios = axios
    // ^ ^ ^ this will allow you to use this.$axios (for Vue Options API form)
    //       so you won't necessarily have to import axios in each vue file

    app.config.globalProperties.$api = api
    // ^ ^ ^ this will allow you to use this.$api (for Vue Options API form)
    //       so you can easily perform requests against your app's API
})

export {api}
