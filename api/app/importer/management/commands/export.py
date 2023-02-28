from django.core.management.base import BaseCommand
from friendlylog import colored_logger as log
import yaml
from django.apps import apps


def export_data(app_name, model_name, filter_field, filter_value, output_file):
    model = apps.get_model(app_label=app_name, model_name=model_name)
    if filter_field and filter_value:
        queryset = model.objects.filter(**{filter_field: filter_value})
    else:
        queryset = model.objects.all()

    data = []
    for obj in queryset:
        fields = {}
        for field in obj._meta.fields:
            fields[field.name] = str(field.value_from_object(obj))
        data.append(fields)

    with open(output_file, 'w') as file:
        yaml.dump(data, file)

    log.info(f'Data exported successfully to {output_file}')


class Command(BaseCommand):

    help = 'Export data from a Django model to a YAML file'

    def add_arguments(self, parser):
        parser.add_argument('app', type=str, help='The app containing the '
                                                  'model to export')
        parser.add_argument('model', type=str, help='The name of the model '
                                                    'to export')
        parser.add_argument('--filter', type=str, help='The filter field '
                                                       'to apply')
        parser.add_argument('--filter-value', type=str, help='The value of '
                                                             'the filter '
                                                             'field')
        parser.add_argument('--output-file', type=str, default='data.yaml',
                            help='The name of the output file')

    def handle(self, *args, **options):
        app_name = options.get('app')
        model_name = options.get('model')
        output_file = options.get('output_file')
        filter_field = options.get('filter')
        filter_value = options.get('filter_value')
        if filter_value and filter_field:
            export_data(app_name, model_name, filter_field, filter_value,
                        output_file)
        else:
            export_data(app_name, model_name, None, None, output_file)