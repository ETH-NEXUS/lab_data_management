import {Notify} from 'quasar'

const handleError = (err: unknown, notify = true) => {
  if (notify) {
    Notify.create({
      message: String(err),
      multiLine: true,
      caption: 'Error',
      icon: 'error',
      color: 'negative',
    })
  }
  console.error(err)
}

export {handleError}
