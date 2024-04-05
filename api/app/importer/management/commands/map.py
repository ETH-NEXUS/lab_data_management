import os
import traceback
from os.path import join
import yaml
from django.core.management.base import BaseCommand
from importer.mappers import EchoMapper, M1000Mapper, MicroscopeMapper, DatMapper
from core.models import Experiment
from importer.helper import message
from core.config import Config
from helpers.logger import logger


def has_csv_files(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".csv"):
                return True
    return False


def die(message):
    raise Exception(message)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "machine",
            type=str,
            choices=("echo", "m1000", "microscope", "dat"),
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
            "--measurement_name", "-n", help="You need to provide a measurement name "
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
                # if in the folder which was provided as 'path' arguments there are no .csv files we use xml_blob, otherwise the file_blob
                if has_csv_files(path):
                    file_blob = Config.current.importer.echo.default.file_blob
                else:
                    file_blob = Config.current.importer.echo.default.xml_blob
                mapper.run(
                    join(path, file_blob),
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
                if not options.get("experiment_name", None):
                    die(
                        "No experiment name provided. If you would like to add missing "
                        "plates, you need to provide the experiment name."
                    )
                measurement_name = options.get("measurement_name", None)
                mapper = M1000Mapper()
                mapper.run(
                    join(path, Config.current.importer.m1000.default.file_blob),
                    debug=options.get("debug", False),
                    measurement_name=measurement_name,
                    experiment_name=options.get("experiment_name", None),
                    room_name=options.get("room_name", None),
                )

            except Exception as ex:
                message(f"Error: {ex}", "error", options.get("room_name", None))
                traceback.print_exc()
        elif options.get("machine") == "microscope":
            try:
                if not options.get("experiment_name", None):
                    die(
                        "No experiment name provided. If you would like to add missing "
                        "plates, you need to provide the experiment name."
                    )
                print(options.get("machine"))
                mapper = MicroscopeMapper()
                mapper.run(
                    join(path, Config.current.importer.microscope.default.file_blob),
                    debug=options.get("debug", False),
                    experiment_name=options.get("experiment_name", None),
                    room_name=options.get("room_name", None),
                )

            except Exception as ex:
                message(f"Error: {ex}", "error", options.get("room_name", None))
                traceback.print_exc()
        elif options.get("machine") == "dat":
            try:
                if not options.get("experiment_name", None):
                    die(
                        "No experiment name provided. If you would like to add missing "
                        "plates, you need to provide the experiment name."
                    )

                mapper = DatMapper()
                mapper.run(
                    join(path, Config.current.importer.dat.default.file_blob),
                    debug=options.get("debug", False),
                    experiment_name=options.get("experiment_name", None),
                    room_name=options.get("room_name", None),
                )

            except Exception as ex:
                message(f"Error: {ex}", "error", options.get("room_name", None))
                traceback.print_exc()
