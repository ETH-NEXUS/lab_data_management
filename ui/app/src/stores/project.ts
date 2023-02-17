import {defineStore} from 'pinia'
import {Project, Experiment, Plate} from 'src/components/models'
import {api} from 'src/boot/axios'
import {ref} from 'vue'

export const useProjectStore = defineStore('project', () => {
  const projects = ref<Array<Project>>([])

  const initialize = async () => {
    const resp_p = await api.get('/api/projects/')
    projects.value = resp_p.data.results
  }

  const add = async (projectName: string) => {
    const resp = await api.post('/api/projects/', {
      name: projectName,
    })
    const project = resp.data
    projects.value.push(project)
    return project
  }

  const addExperiment = async (project: Project, experimentName: string) => {
    const resp = await api.post('/api/experiments/', {
      name: experimentName,
      project: project.id,
    })
    const experiment = resp.data
    project.experiments.push(experiment)
    return experiment
  }

  const addPlate = async (experiment: Experiment, barcode: string) => {
    const resp = await api.post('/api/plates/', {
      barcode: barcode,
      experiment: experiment.id,
    })
    const plate = resp.data
    experiment.plates.push(plate)
    return plate
  }

  const generateBarcodes = async (
    experimentId: number,
    prefix: string,
    numberOfPlates: number,
    sides: Array<string>
  ) => {
    await api.post('/api/experiments/barcodes/', {
      number_of_plates: numberOfPlates,
      experiment_id: experimentId,
      sides: sides,
      prefix: prefix,
    })
  }

  const updateBarcode = async (id: number, prefix: string, numberOfPlates: number, sides: Array<string>) => {
    console.log(id, prefix, numberOfPlates, sides)
    await api.patch(`/api/barcodespecifications/${id}/`, {
      number_of_plates: numberOfPlates,
      id: id,
      sides: sides,
      prefix: prefix,
    })
  }

  const deleteBarcode = async (id: number) => {
    await api.delete(`/api/barcodespecifications/${id}/`)
  }

  return {projects, initialize, add, addExperiment, addPlate, generateBarcodes, updateBarcode, deleteBarcode}
})
