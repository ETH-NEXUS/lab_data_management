import {boot} from 'quasar/wrappers'
import VueExcelEditor from 'vue3-excel-editor'

export default boot(({app}) => {
  app.use(VueExcelEditor)
})
