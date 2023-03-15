import {defineStore} from 'pinia'
import {Experiment, PlateDimension, Project, ExperimentPayload} from 'src/components/models'
import {api} from 'src/boot/axios'
import {ref} from 'vue'

export const useProjectStore = defineStore('project', () => {
  const projects = ref<Array<Project>>([])
  const plateDimensions = ref<Array<PlateDimension>>([])
  const experiments = ref<Array<Experiment>>([])

  const initialize = async () => {
    const resp_p = await api.get('/api/projects/')
    projects.value = resp_p.data.results

    const resp_e = await api.get('/api/experiments/')
    experiments.value = resp_e.data.results

    const res_d = await api.get('/api/platedimensions/')
    plateDimensions.value = res_d.data.results
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

  const addPlatesToExperiment = async (
    experimentId: number,
    barcodeSpecificationsId: number,
    dimension: number
  ) => {
    alert(
      'experimentId: ' +
        experimentId +
        ' barcodeSpecificationsId: ' +
        barcodeSpecificationsId +
        ' dimension: ' +
        dimension +
        ''
    )
    await api.post('/api/experiments/bulk_add_plates/', {
      experiment_id: experimentId,
      barcode_specification_id: barcodeSpecificationsId,
      plate_dimension_id: dimension,
    })
  }

  const updateExperiment = async (experimentId: number, payload: ExperimentPayload) => {
    await api.patch(`/api/experiments/${experimentId}`, payload)
  }

  return {
    projects,
    initialize,
    add,
    addExperiment,
    updateExperiment,
    addPlate,
    generateBarcodes,
    updateBarcode,
    deleteBarcode,
    plateDimensions,
    addPlatesToExperiment,
    experiments,
  }
})
