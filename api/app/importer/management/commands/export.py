from django.core.management.base import BaseCommand, CommandError
from django.core import serializers
from django.apps import apps
import os
from os import makedirs
from helpers.logger import logger


def export_data(app_name, model_name, filters=None):
    """
    The objects of the model as YAML. Each filter is "field=value", e.g.
    ["barcode__startswith=Drug08", "library__name=LLD_24000"].
    """
    model = apps.get_model(app_label=app_name, model_name=model_name)
    queryset = model.objects.all()
    for condition in filters or []:
        if "=" not in condition:
            raise CommandError(f"The filter '{condition}' is not field=value.")
        field, value = condition.split("=", 1)
        queryset = queryset.filter(**{field: value})

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
            "--filter",
            "-f",
            action="append",
            help="A filter field=value, e.g. barcode__startswith=Drug08; can be repeated",
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
        data = export_data(app_name, model_name, options.get("filter"))

        makedirs(os.path.split(output_file)[0], exist_ok=True)
        with open(output_file, "a" if options.get("append") else "w") as file:
            file.write(data)

        logger.info(f"Data exported successfully to {output_file}")
