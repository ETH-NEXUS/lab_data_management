import traceback
from importer.mapping import EchoMapping, MeasurementMapping
from django.core.management.base import BaseCommand
from friendlylog import colored_logger as log
import yaml


# add plates to the experiment, put the echo files to /data/echo and run the
# following command:

# ./manage.py map -t echo -p /data/echo -m /data/echo/echo-headers.yml -e exp

# you can try it out with the LLD compound library. Make sure to add plates
# with the barcode prefix BAF210901 and the correct dimensions to the
# experiment first. This can be done in the UI.


# ./manage.py map -t measurement -p /data/M1000  -e exp

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--type', '-t', type=str, required=True,
                            help='Type of mapping, e. g. echo')
        parser.add_argument('--path', '-p', type=str, required=True,
                            help='Path to the directory containing the '
                                 'mapping files')
        parser.add_argument('--experiment', '-e', type=str, required=True)
        parser.add_argument('--mapping-file', '-m', type=str,
                            help='A yml file with the column headers, otherwise default '
                                 'headers are used')

    def handle(self, *args, **options):
        path = options.get('path')
        if options.get('type') == 'echo':
            headers = EchoMapping.DEFAULT_COLUMN_HEADERS
            headers_file = options.get('headers_file', None)
            if headers_file:
                try:
                    with open(options.get('headers_file'), 'r') as file:
                        headers = yaml.safe_load(file)
                except FileNotFoundError:
                    log.error(
                        f"The headers file '{headers_file}' could not be found.")
                    return
                except yaml.YAMLError as e:
                    log.error(
                        f"Error parsing the YAML file '{headers_file}': {e}")
                    return

            try:
                data = EchoMapping.get_csv_echo_files(path, headers)
                for mapping in data:
                    EchoMapping.parse_plate_data(mapping['data'],
                                                 mapping['file_path'], headers,
                                                 options.get('experiment'))

            except Exception as ex:
                log.error(ex)
                traceback.print_exc()

        elif options.get('type') == 'measurement':
            data = MeasurementMapping.get_measurement_files(path)
            for item in data:
                try:
                    MeasurementMapping.parse_measurement_data(item['measurement_data'], item['barcode'])
                except Exception as ex:
                    log.error(ex)
                    traceback.print_exc()
