import {defineStore} from 'pinia'

import {
  Experiment,
  PlateDimension,
  Project,
  ExperimentPayload,
  ProjectPayload,
  harvestProject,
} from 'src/components/models'

import {api} from 'src/boot/axios'
import {ref} from 'vue'

export const useProjectStore = defineStore('project', () => {
  const projects = ref<Array<Project>>([])
  const plateDimensions = ref<Array<PlateDimension>>([])
  const experiments = ref<Array<Experiment>>([])
  const harvestProjects = ref<harvestProject[]>([])

  const initialize = async () => {
    const resp_p = await api.get('/api/projects/')
    projects.value = resp_p.data.results

    const resp_e = await api.get('/api/experiments/')
    experiments.value = resp_e.data.results

    const res_d = await api.get('/api/platedimensions/')
    plateDimensions.value = res_d.data.results

    const resp_harvest = await api.get('/api/harvest/harvest_projects/')
    const harvestDataCopy = JSON.parse(JSON.stringify(resp_harvest.data.projects))
    for (let i = 0; i < harvestDataCopy.length; i++) {
      harvestDataCopy[i].value = harvestDataCopy[i].name
      harvestDataCopy[i].label = harvestDataCopy[i].name
    }
    harvestProjects.value = harvestDataCopy.filter((item: harvestProject) => item.is_active)
    console.log(harvestProjects.value)
  }

  const add = async (payload: {projectName: string; harvest_id: number | null}) => {
    const resp = await api.post('/api/projects/', payload)
    const project = resp.data
    projects.value.push(project)
    return project
  }
  const updateHarvestInfo = async (projectId: number) => {
    const resp = await api.get(`/api/harvest/update_harvest_info/${projectId}/`)
    await initialize()
    return resp
  }

  const updateProject = async (projectId: number, payload: ProjectPayload) => {
    await api.patch(`/api/projects/${projectId}`, payload)
    await initialize()
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
    initialize()
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
    updateProject,
    harvestProjects,
    updateHarvestInfo,
  }
})
