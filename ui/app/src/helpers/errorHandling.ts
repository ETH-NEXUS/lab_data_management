import {Notify} from 'quasar'
import {AxiosError} from 'axios'

export const formatKV = (obj: object) => {
  if (!obj) {
    return 'No details available'
  }
  let ret = ''
  for (const [key, value] of Object.entries(obj)) {
    if (ret !== '') {
      ret += ' / '
    }
    ret += `${key}: ${value}`
  }
  return ret
}

export const handleError = (err: AxiosError | string | unknown, notify = true) => {
  if (notify) {
    let caption = 'error.no_details_available'
    if (err instanceof Error) {
      if ((err as AxiosError).response?.headers['content-type'].startswith('text/html')) {
        caption = (err as AxiosError).response?.statusText || 'error.no_details_available'
      } else {
        caption = (err as AxiosError).response?.data.detail || formatKV((err as AxiosError).response?.data)
      }
    }
    Notify.create({
      message: String(err),
      multiLine: true,
      caption: caption,
      icon: 'error',
      color: 'negative',
      timeout: 0,
      closeBtn: true,
    })
  }
  console.error(err)
}

export const success = (msg: string, caption = '') => {
  Notify.create({
    message: msg,
    multiLine: true,
    caption: caption,
    icon: 'info',
    color: 'positive',
  })
}
