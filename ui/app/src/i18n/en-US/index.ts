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
    notebook: 'Notebook',
    about: 'About',
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
    chain: 'Chain',
    well_type: 'Type',
    template_name: 'Template Name',
    status: 'Status',
  },
  message: {
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

    custom_name: 'Custom name',
  },
  action: {
    add_plates: 'Add paltes',
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
  },
  error: {
    select_plate_dimension: 'Please select a plate dimension.',
    no_details_available: 'No details available.',
  },
  unit: {
    amount: 'µl',
  },
  info: {
    no_description: 'No description provided',
    applying_in_progress: 'Applying template in progress. Hang on...',
  },
}
