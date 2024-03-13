// This is just an example,
// so you can safely delete all default props below

// title: field === 'name' ? t('project.edit_project_name') : t('project.edit_project_description'),
// message: field === 'name' ? t('project.project_name') : t('project.project_description'),

export default {
  project: {
    project_description: 'Project description',
    number_of_experiments: 'Number of experiments',
    project_name: 'Project name',
    created_at: 'Created at',
    no_description: 'No description provided',
    edit_project_name: 'Edit project name',
    edit_project_description: 'Edit project description',
    harvest_notes: 'Harvest notes',
  },

  experiment: {
    description: 'Description',
    number_plates: 'Number of plates',
    created_at: 'Created at',
    barcode_sets: 'Barcode sets',
    delete_specifications: 'Delete specifications',
    edit_specifications: 'Edit specifications',
    add_plates: 'Add plates to the experiment using these specifications',
    add_plates_to_experiment: 'Add plates to the experiment',
    choose_dimensions: 'Choose the dimensions of the plate',
    plates_added:
      'You have added the plates with these barcode specifications to this\n' +
      '                        experiment',
  },
  failed: 'Action failed',
  success: 'Action was successful',
  label: {
    username: 'Username',
    password: 'Password',
    login: 'Login',
    logout: 'Logout',
    compound_libraries: 'Compound Libraries',
    projects: 'Projects',
    filter: 'Filter',
    name: 'Name',
    abbrev: 'Abbrev',
    value: 'Value',
    unit: 'Unit',
    identifier: 'Identifier',
    structure: 'Structure',
    amount: 'Amount',
    initial_amount: 'Initial Amount',
    plate_dimension: 'Plate Dimension',
    position: 'Position',
    compound: 'Compound',
    ok: 'Ok',
    add: 'Add',
    cancel: 'Cancel',
    plate: 'Plate',
    well: 'Well',
    timestamp: 'Timestamp',
    target_well: 'Target well',
    source_well: 'Source well',
    target_plate: 'Target plate',
    map: 'Map',
    copy: 'Copy',
    delimiter: 'Delimiter',
    quotechar: 'Quotechar',
    from_column: 'From column',
    to_column: 'To column',
    amount_column: 'Amount column',
    measurement_feature: 'Measurement feature',
    measurement_value: 'Measurement value',
    templates: 'Templates',
    template_plate: 'Template plate',
    well_content: 'Well content',
    hr_position: 'Position',
    index_position: 'Index',
    well_type: 'Well type',
    show_heatmap: 'Show heatmap',
    select_measurement: 'Select measurement',
    select_color_palette: 'Select color palette',
    compounds: 'Compounds',
    measurements: 'Measurements',
    z_prime: "z'",
    ssmd: 'ssmd',
    notebook: 'Notebook',
    about: 'About',
    help: 'Help',
    smaller_map_view: 'Smaller map view',
    label: 'Label',
    plot_view: 'Plot view',
    show_type_as_square: 'Show wells as squares',
  },
  hint: {
    username: 'Enter your ETH username',
    password: 'Enter your ETH password',
    compound_to_add: 'Select the compound to add',
    well_to_transfer_from: 'Select the well to transfer the compound from',
    amount_to_transfer: 'Enter the amount to transfer',
    target_plate: 'Enter the barcode of the target plate you want to map this plate to',
    measurement_feature: 'Select the mMeasurement feature',
    measurement_value: 'Enter the measurement value',
    template_plate: 'Select a template plate',
  },
  title: {
    app: 'Lab Data Management',
    compound: 'Compound',
    compounds: 'Compounds',
    well: 'Well',
    measurements: 'Measurements',
    project_name: 'Project Name',
    experiment_name: 'Experiment Name',
    edit_experiment_name: 'Edit Experiment Name',
    edit_experiment_description: 'Edit Experiment Description',
    plate_barcode: 'Plate Barcode',
    withdrawals: 'Withdrawals',
    donors: 'Donors',
    amount: 'Amount',
    amount_dmso: 'Amount &  DMSO',
    chain: 'Chain',
    well_type: 'Type',
    template_name: 'Template Name',
    status: 'Status',
    sure_delete: 'Are you sure you want to delete this item?',
  },
  message: {
    select_label: 'Select label',
    select_timestamp: 'Select measurement time',
    successfully_logged_in: 'Successfully logged in',
    logged_out: 'You are logged out - see you soon again',
    project_name: 'Please enter the project name:',
    experiment_name: 'Please enter the experiment name:',
    experiment_description: 'Please enter the experiment description:',
    plate_barcode: 'Please enter (or scan) the plates barcode:',
    plate_has_no_dimension: 'The plate has no dimension assigned yet. Please specify one.',
    no_well_information: 'No existing well information do you want to create it?',
    no_compounds_found: 'No compounds found',
    no_plates_found: 'No plates found',
    no_wells_found: 'No wells found',
    successfully_mapped_plate: 'Successfully mapped plate',
    successfully_copied_plate: 'Successfully copied plate',
    no_dimension: 'No dimension',
    template_name: 'Please enter the template name:',
    select_dimension: 'Please select a dimension',
    project_name_harvest: 'Select project name from Harvest',
    update_harvest: 'Update Harvest information',
    harvest_info_updated: 'Harvest information updated',
    project_name_required: 'Project name is required. Please enter something.',

    custom_name: 'Custom name',
  },
  action: {
    update_harvest_projects: 'Update Harvest projects',
    calculate_measurement: 'New measurement',
    add_plates: 'Add paltes',
    show_results: 'Show measurement results',
    add_barcode_specification: 'Add barcode specification',
    new_project: 'New project',
    project_properties: 'Project properties',
    new_experiment: 'New Experiment',
    submit: 'Submit',
    experiment_properties: 'Experiment properties',
    generate_barcodes: 'Generate barcodes',
    enter_prefix: 'Enter a prefix *',
    enter_number_of_plates: 'Enter the number of plates *',
    select_sides: 'Select sides',
    validation_number_of_plates: 'Please enter a number',
    validation_positive_number_of_plates: 'Please enter a positive number',
    validation_prefix: 'Please type something',
    new_plate: 'New Plate',
    create_well: 'Create Well',
    add_compound: 'Add Compound',
    add_measurement: 'Add Measurement',
    map_plate: 'Map Plate',
    copy_plate: 'Copy Plate',
    compoundlib_properties: 'Compound Library Properties',
    new_category: 'New Category',
    category_properties: 'Category Properties',
    new_template: 'New Template',
    apply_template: 'Apply Template',
    apply: 'Apply',
    delete_specification: 'Delete specification',
    edit_specification: 'Edit specification',
    download_csv: 'Download CSV',
    delete_item: 'Delete item',
    download_file: 'Download file',
    upload_file: 'Upload file to the directory',
    view_file_content: 'View file content',
    delete_well: 'Delete well',
    generate_report: 'Generate report',
    download_csv_data: 'Download CSV data',
  },
  error: {
    select_plate_dimension: 'Please select a plate dimension.',
    no_details_available: 'No details available.',
    cannot_login: 'Cannot login.',
  },
  unit: {
    mikro: 'µl',
    nL: 'nL',
  },
  info: {
    downloading_results: 'Downloading results...',
    no_description: 'No description provided',
    applying_in_progress: 'Applying template in progress. Hang on...',
    calculation_in_progress: 'Calculation in progress. Hang on...',
    running_in_progress: 'Running in progress. Hang on...',
    generation_in_progress: 'Generation in progress. Hang on...',
  },
  management: {
    selected_directories: 'Selected directories/files',
    no_dirs_selected: 'No directories/files selected',
    management: 'Management',
    echo: 'Echo/M1000/Micro',
    import: 'Import',
    export: 'Export',
    import_control_plate: 'Import control plate',
    import_library_plate: 'Import library plate',
    import_sdf_library: 'Import SDF library',
  },
  well: {
    current_amount: 'Current amount',
    current_dmso: 'Current  DMSO',
  },
}
