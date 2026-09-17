import os
from os.path import join
import yaml
from django.core.management.base import BaseCommand, CommandError
from importer.mappers import EchoMapper, M1000Mapper, MicroscopeMapper
from core.models import Experiment
from importer.helper import message
from importer.command_output import error_text
from importer.config import Config
from helpers.logger import logger


def has_csv_files(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".csv"):
                return True
    return False


def die(message):
    raise CommandError(message)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "machine",
            type=str,
            choices=("echo", "m1000", "microscope", "C10-reader", "C10-imager"),
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

    @staticmethod
    def read_echo_columns(path: str) -> dict:
        """
        The Echo column names from a yml file, with the same keys as `columns`
        of echo in ldm.yaml, e.g. {"source_well": "Source Well", "DMSO": "% DMSO", ...}.
        """
        try:
            with open(path, "r") as file:
                columns = yaml.safe_load(file)
        except FileNotFoundError:
            raise CommandError(f"The column file '{path}' could not be found.")
        except yaml.YAMLError as error:
            raise CommandError(f"Error parsing the YAML file '{path}': {error}")

        if not isinstance(columns, dict):
            raise CommandError(
                f"The column file '{path}' has no column names, "
                "e.g. 'source_well: Source Well'."
            )
        missing_keys = [key for key in EchoMapper.DEFAULT_COLUMNS if key not in columns]
        if missing_keys:
            raise CommandError(
                f"The column file '{path}' is missing these keys: "
                f"{', '.join(missing_keys)}."
            )
        return columns

    @staticmethod
    def show_error(error: Exception, options: dict) -> None:
        """Shows the error on the management page; an unexpected one also goes to the log."""
        message(error_text(error), "error", options.get("room_name"))
        if not isinstance(error, CommandError):
            logger.exception(f"Command map {options.get('machine')} failed")

    def handle(self, *args, **options):

        path = options.get("path")
        if not os.path.exists(path):
            raise CommandError(f"The folder {path} does not exist.")
        if not os.path.isdir(path):
            raise CommandError(
                f"{path} is a file. Please choose the folder that contains it."
            )
        if options.get("experiment_name", None):
            experiment = Experiment.objects.filter(
                name=options.get("experiment_name")
            ).first()
            if not experiment:
                raise CommandError(
                    f"No experiment with name '{options.get('experiment_name')}' found in the database."
                )

        if options.get("machine") == "echo":
            headers = EchoMapper.DEFAULT_COLUMNS
            if options.get("mapping_file"):
                headers = self.read_echo_columns(options.get("mapping_file"))
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
            except Exception as error:
                self.show_error(error, options)

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

            except Exception as error:
                self.show_error(error, options)
        elif options.get("machine") in ["microscope", "C10-imager", "C10-reader"]:
            try:
                if not options.get("experiment_name", None):
                    die(
                        "No experiment name provided. If you would like to add missing "
                        "plates, you need to provide the experiment name."
                    )
                mapper = MicroscopeMapper()

                if path.endswith(".txt"):
                    _blob = Config.current.importer.microscope.default.file_blob
                else:
                    _blob = Config.current.importer.microscope.default.txt_blob
                    logger.info(f"Using blob: {_blob}")
                mapper.run(
                    join(path, _blob),
                    debug=options.get("debug", False),
                    experiment_name=options.get("experiment_name", None),
                    room_name=options.get("room_name", None),
                    measurement_name=options.get("measurement_name", None),
                )

            except Exception as error:
                self.show_error(error, options)
