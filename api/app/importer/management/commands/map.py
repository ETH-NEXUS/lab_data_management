import traceback
from os.path import join
import yaml
from django.core.management.base import BaseCommand
from importer.mappers import EchoMapper, M1000Mapper
from core.models import Experiment
from importer.helper import message
from core.config import Config


def die(message):
    raise Exception(message)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "machine",
            type=str,
            choices=("echo", "m1000"),
            help="Machine to map from",
        )
        parser.add_argument(
            "--path",
            "-p",
            type=str,
            required=True,
            help="Path to the directory containing the " "mapping files",
        )
        parser.add_argument(
            "--mapping-file",
            "-m",
            type=str,
            help="A yml file with the column headers, "
            "otherwise default headers are used",
        )
        parser.add_argument(
            "--debug",
            "-d",
            action="store_true",
            help="Enable debug mode",
        )

        parser.add_argument(
            "--evaluate",
            "-e",
            help="A formula which has to be applied to the "
            "measurement "
            "values, e. g. --evaluate 'Acceptor/Donor'",
        )
        parser.add_argument(
            "--measurement_name",
            "-n",
            help="If an evaluation formula is provided, "
            "you need to provide a measurement name "
            "as well'",
        )

        parser.add_argument(
            "--create_missing_plates",
            "-c",
            action="store_true",
            help="Create missing plates by measurement mapping",
        )

        parser.add_argument(
            "--room_name",
            "-r",
            help="Unique room name for long polling.",
        )
        parser.add_argument(
            "--experiment_name",
            "-x",
            help="If you would like to create missing plates by measurement "
            "mapping, you need to"
            "provide the experiment name.",
        )

    def handle(self, *args, **options):
        path = options.get("path")
        if options.get("experiment_name", None):
            experiment = Experiment.objects.filter(
                name=options.get("experiment_name")
            ).first()
            if not experiment:
                message(
                    f"No experiment with name '{options.get('experiment_name')}' found in the database.",
                    "error",
                    options.get("room_name", None),
                )
                raise ValueError(
                    f"No experiment with name '{options.get('experiment_name')}' found in the database."
                )

        if options.get("machine") == "echo":
            headers = EchoMapper.DEFAULT_COLUMNS
            headers_file = options.get("headers_file", None)

            if headers_file:
                try:
                    with open(options.get("headers_file"), "r") as file:
                        headers = yaml.safe_load(file)
                except FileNotFoundError:
                    message(
                        f"The headers file '{headers_file}' could not be found.",
                        "error",
                        options.get("room_name", None),
                    )
                    return
                except yaml.YAMLError as e:
                    message(
                        f"Error parsing the YAML file '{headers_file}': {e}",
                        "error",
                        options.get("room_name", None),
                    )
                    return
            try:
                mapper = EchoMapper()
                mapper.run(
                    # join(path, "**", "*[_-][Tt]ransfer[_-]*.csv"),
                    join(path, Config.current.importer.echo.default.file_blob),
                    headers=headers,
                    debug=options.get("debug", False),
                    room_name=options.get("room_name", None),
                    experiment_name=options.get("experiment_name", None),
                )
            except Exception as ex:
                message(f"Error: {ex}", "error", options.get("room_name", None))
                traceback.print_exc()

        elif options.get("machine") == "m1000":
            try:
                if options.get("create_missing_plates", False):
                    if not options.get("experiment_name", None):
                        die(
                            "No experiment name provided. If you would like to add missing "
                            "plates, you need to provide the experiment name."
                        )

                evaluation = options.get("evaluate", None)
                measurement_name = options.get("measurement_name", None)
                if evaluation and not measurement_name:
                    raise ValueError("No measurement name provided")
                mapper = M1000Mapper()
                mapper.run(
                    join(path, Config.current.importer.m1000.default.file_blob),
                    debug=options.get("debug", False),
                    evaluation=evaluation,
                    measurement_name=measurement_name,
                    create_missing_plates=options.get("create_missing_plates", False),
                    experiment_name=options.get("experiment_name", None),
                    room_name=options.get("room_name", None),
                )

            except Exception as ex:
                message(f"Error: {ex}", "error", options.get("room_name", None))
                traceback.print_exc()
