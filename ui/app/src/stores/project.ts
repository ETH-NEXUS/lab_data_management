import {defineStore} from 'pinia'

import {
  Experiment,
  PlateDimension,
  Project,
  ExperimentPayload,
  ProjectPayload,
  CalculatorPayload,
  PlateInfo,
  Plate,
} from 'src/components/models'

import {api} from 'src/boot/axios'
import {ref} from 'vue'

export const useProjectStore = defineStore('project', () => {
  const projects = ref<Array<Project>>([])
  const plateDimensions = ref<Array<PlateDimension>>([])
  const experiments = ref<Array<Experiment>>([])
  const outputNotebooks = ref<Array<string>>([])
  const inputNotebookPath = ref<string | null>(null)
  const prefilledPlateInfo = ref<PlateInfo[]>([])
  const controlPlates = ref<Array<Plate>>([])

  const initialize = async () => {
    const resp_p = await api.get('/api/projects/')
    projects.value = resp_p.data.results

    const resp_e = await api.get('/api/experiments/')
    experiments.value = resp_e.data.results

    const res_d = await api.get('/api/platedimensions/')
    plateDimensions.value = res_d.data.results
  }

  const add = async (payload: {name: string; harvest_id: number | null}) => {
    const resp = await api.post('/api/projects/', payload)
    const project = resp.data
    projects.value.push(project)
    return project
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

  const getExperimentPlates = async (experimentId: number) => {
    const resp = await api.get(`/api/plates/?experiment=${experimentId}`)
    return resp.data.results
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

  //   const add_new_measurement = (expression: string, newLabel:string) => {
  //   projectStore.calculateNewMeasurement(props.plate.id, expression, newLabel)
  // }
  const addNewMeasurement = async (
    plateId: number | null,
    expression: string,
    newLabel: string,
    usedLabels: string[],
    experimentId: null | number | undefined = null
  ) => {
    const payload: CalculatorPayload = {
      expression: expression,
      new_label: newLabel,
      used_labels: usedLabels,
      separate_time_series_points: false,
    }
    if (plateId) {
      payload.plate_id = plateId
    } else if (experimentId) {
      payload.experiment_id = experimentId
    }
    if (usedLabels.some(label => label.includes('-->'))) {
      payload.separate_time_series_points = true
    }

    const res = await api.post(`/api/plates/${plateId}/add_new_measurement/`, payload)
    if (res.status === 200) {
      await initialize()
    }
    return res.status
  }

  const generateReport = async (experiment: string, label: string, notebook_path: string) => {
    const res = await api.post('/api/generate_pdf_report/', {
      label: label,
      experiment: experiment,
      notebook_path: notebook_path,
    })
    await getNotebookOutputFiles(experiment)
    return res.data
  }

  const getNotebookOutputFiles = async (experiment: string) => {
    const res = await api.post('/api/list_files/', {experiment: experiment})
    console.log('notebooks', res)
    outputNotebooks.value = res.data.notebooks
  }

  const downloadPDFReport = async (path: string) => {
    const res = await api.post('/api/download_pdf_report/', {path: path}, {responseType: 'blob'})

    console.log('download', res)
    const link = document.createElement('a')
    link.href = window.URL.createObjectURL(new Blob([res.data]))
    link.setAttribute('download', path.split('/').pop() as string)
    document.body.appendChild(link)
    link.click()
  }

  const downloadCSVData = async (experiment: string, label: string) => {
    const res = await api.post('/api/download_csv_data/', {
      label: label,
      experiment: experiment,
    })
    const link = document.createElement('a')
    link.href = window.URL.createObjectURL(new Blob([res.data]))
    link.setAttribute('download', `${experiment}_${label}.csv`)
    document.body.appendChild(link)
    link.click()
  }

  const getPrefilledPlateInfo = async (experimentId: number) => {
    const res = await api.get(`/api/prefillPlateInfo/?experiment_id=${experimentId}`)
    prefilledPlateInfo.value = res.data.plate_info
  }

  const savePlateInfo = async (experimentId: number, plateInfo: PlateInfo[]) => {
    const res = await api.post('/api/save_plate_info/', {
      experiment_id: experimentId,
      plate_info: plateInfo,
    })
    if (res.status === 200) {
      await getPrefilledPlateInfo(experimentId)
      return 'success'
    }
  }

  const getControlPlates = async () => {
    const res = await api.get('/api/plates/?is_control_plate=true')
    controlPlates.value = res.data.results
  }

  return {
    projects,
    plateDimensions,
    experiments,
    outputNotebooks,
    inputNotebookPath,
    prefilledPlateInfo,
    controlPlates,
    initialize,
    add,
    addExperiment,
    updateExperiment,
    addPlate,
    generateBarcodes,
    updateBarcode,
    deleteBarcode,
    updateProject,
    addNewMeasurement,
    getExperimentPlates,
    generateReport,
    getNotebookOutputFiles,
    downloadPDFReport,
    downloadCSVData,
    getPrefilledPlateInfo,
    savePlateInfo,
    getControlPlates,
  }
})
