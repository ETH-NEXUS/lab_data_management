import traceback
from os.path import join

import yaml
from django.core.management.base import BaseCommand
from friendlylog import colored_logger as log

from importer.mappers import EchoMapper, M1000Mapper


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

    def handle(self, *args, **options):
        path = options.get("path")

        if options.get("machine") == "echo":
            headers = EchoMapper.DEFAULT_COLUMNS
            headers_file = options.get("headers_file", None)
            if headers_file:
                try:
                    with open(options.get("headers_file"), "r") as file:
                        headers = yaml.safe_load(file)
                except FileNotFoundError:
                    log.error(f"The headers file '{headers_file}' could not be found.")
                    return
                except yaml.YAMLError as e:
                    log.error(f"Error parsing the YAML file '{headers_file}': {e}")
                    return
            try:
                mapper = EchoMapper()
                mapper.run(
                    join(path, "**", "*-transfer-*.csv"),
                    headers=headers,
                    debug=options.get("debug", False),
                )
            except Exception as ex:
                log.error(ex)
                traceback.print_exc()

        elif options.get("machine") == "m1000":
            try:
                evaluation = options.get("evaluate", None)
                measurement_name = options.get("measurement_name", None)
                if evaluation and not measurement_name:
                    raise ValueError("No measurement name provided")
                mapper = M1000Mapper()
                mapper.run(
                    join(path, "*.asc"),
                    debug=options.get("debug", False),
                    evaluation=evaluation,
                    measurement_name=measurement_name,
                )

            except Exception as ex:
                log.error(ex)
                traceback.print_exc()
