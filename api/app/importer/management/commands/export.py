from django.core.management.base import BaseCommand
from django.core import serializers
from django.apps import apps
import os
from os import makedirs


def export_data(app_name, model_name, _filter=None):
    model = apps.get_model(app_label=app_name, model_name=model_name)
    if _filter:
        # TODO: Introduce a flexible filter mechanism
        print(_filter)
        queryset = model.objects.filter(eval(_filter))
    else:
        queryset = model.objects.all()

    return serializers.serialize("yaml", queryset)


class Command(BaseCommand):
    help = "Export data from a Django model to a YAML file"

    def add_arguments(self, parser):
        parser.add_argument(
            "app", type=str, help="The app containing the " "model to export"
        )
        parser.add_argument(
            "model", type=str, help="The name of the " "model " "to export"
        )
        parser.add_argument(
            "--filter", "-f", type=str, help="The filter " "field " "to apply"
        )
        parser.add_argument(
            "--append",
            "-a",
            action="store_true",
            help="Append the result to the existing file",
        )

        parser.add_argument(
            "--output-file",
            "-o",
            type=str,
            default="test_data/data.yaml",
            help="The name of the output file",
        )

    def handle(self, *args, **options):
        app_name = options.get("app")
        model_name = options.get("model")
        output_file = options.get("output_file")
        _filter = options.get("filter")
        if _filter and not _filter.startswith("Q("):
            _filter = f"Q({_filter})"
        data = export_data(app_name, model_name, _filter)

        makedirs(os.path.split(output_file)[0], exist_ok=True)
        with open(output_file, "a" if options.get("append") else "w") as file:
            file.write(data)

        print(f"Data exported successfully to {output_file}")
