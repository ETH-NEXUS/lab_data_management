import os
from os.path import join

import yaml
from django.core.management.base import BaseCommand, CommandError

from core.models import Experiment
from importer.config import Config
from importer.helper import message
from importer.mappers import BaseMapper, EchoMapper, M1000Mapper, MicroscopeMapper

# Without an experiment name the mappers cannot create a missing plate
NO_EXPERIMENT_NAME = (
    "No experiment name provided. If you would like to add missing plates, "
    "you need to provide the experiment name."
)


def first_pattern_with_files(path: str, patterns: tuple[str, ...]) -> str:
    """
    The first file pattern of `patterns` that finds files in the folder, e.g.
    the Echo CSV pattern, or the XML pattern when the folder has no CSV reports.
    The last pattern is used when none of them finds a file, so that the mapper
    can say which files it looked for.
    """
    for pattern in patterns:
        if BaseMapper.get_files(join(path, pattern)):
            return pattern
    return patterns[-1]


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
            "--measurement_name",
            "-n",
            help="The label of the measured values, e.g. 'Lum'. Without it, the "
            "label of the file is used.",
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
    def warn_about_unused_column_file(options: dict) -> None:
        """Only the Echo mapper reads a column file, the other machines ignore it."""
        if options.get("mapping_file"):
            message(
                f"The column file {options['mapping_file']} is only used for echo "
                f"reports, not for {options.get('machine')}.",
                "warning",
                options.get("room_name"),
            )

    @staticmethod
    def warn_about_unused_measurement_name(options: dict) -> None:
        """An Echo report has no measurements, so it has no measurement name."""
        if options.get("measurement_name"):
            message(
                f"The measurement name {options['measurement_name']} is not used "
                "for echo reports, only for measurement files.",
                "warning",
                options.get("room_name"),
            )

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
            self.warn_about_unused_measurement_name(options)
            headers = EchoMapper.DEFAULT_COLUMNS
            if options.get("mapping_file"):
                headers = self.read_echo_columns(options.get("mapping_file"))
            echo = Config.current.importer.echo.default
            pattern = first_pattern_with_files(path, (echo.file_blob, echo.xml_blob))
            EchoMapper().run(
                join(path, pattern),
                headers=headers,
                debug=options.get("debug", False),
                room_name=options.get("room_name"),
                experiment_name=options.get("experiment_name"),
            )

        elif options.get("machine") == "m1000":
            if not options.get("experiment_name"):
                raise CommandError(NO_EXPERIMENT_NAME)
            self.warn_about_unused_column_file(options)
            M1000Mapper().run(
                join(path, Config.current.importer.m1000.default.file_blob),
                debug=options.get("debug", False),
                measurement_name=options.get("measurement_name"),
                experiment_name=options.get("experiment_name"),
                room_name=options.get("room_name"),
            )

        elif options.get("machine") in ["microscope", "C10-imager", "C10-reader"]:
            if not options.get("experiment_name"):
                raise CommandError(NO_EXPERIMENT_NAME)
            self.warn_about_unused_column_file(options)
            # The C10 writes .txt files in the reader mode and .xlsx in the imager
            # mode, so the chosen mode decides which files are read. The old name
            # "microscope" does not say the mode, there both are looked for.
            microscope = Config.current.importer.microscope.default
            if options.get("machine") == "C10-reader":
                pattern = microscope.txt_blob
            elif options.get("machine") == "C10-imager":
                pattern = microscope.file_blob
            else:
                pattern = first_pattern_with_files(
                    path, (microscope.txt_blob, microscope.file_blob)
                )
            MicroscopeMapper().run(
                join(path, pattern),
                debug=options.get("debug", False),
                experiment_name=options.get("experiment_name"),
                room_name=options.get("room_name"),
                measurement_name=options.get("measurement_name"),
            )
