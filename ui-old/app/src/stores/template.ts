import {defineStore} from 'pinia'
import {TemplateCategory, Template} from 'src/components/models'
import {api} from 'src/boot/axios'
import {ref} from 'vue'

export const useTemplateStore = defineStore('template', () => {
  const templateCategories = ref<Array<TemplateCategory>>([])

  const initialize = async () => {
    const resp_p = await api.get('/api/templatecategories/')
    templateCategories.value = resp_p.data.results
  }

  const addCategory = async (templateCategoryName: string) => {
    const resp = await api.post('/api/templatecategories/', {
      name: templateCategoryName,
    })
    const category = resp.data
    templateCategories.value.push(category)
    return category
  }

  const addTemplate = async (category: TemplateCategory, templateName: string) => {
    const respTemplate = await api.post('/api/templates/', {
      name: templateName,
      category: category.id,
    })
    const template = respTemplate.data

    const plate = await addPlate(template, category.name)
    template.plate = plate

    category.templates.push(template)
    return template
  }

  const addPlate = async (template: Template, categoryName: string) => {
    const respPlate = await api.post('/api/plates/', {
      barcode: `__TEMPL__${categoryName}_${template.name}`,
      template: template.id,
    })
    return respPlate.data
  }

  return {templateCategories, initialize, addCategory, addTemplate, addPlate}
})
