import inspect
from importlib import import_module
from pathlib import Path

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from helpers.logger import logger

import_inventory_csv = import_module(
    "inventory.import.services.inventory_csv_importer"
).import_inventory_csv


class Command(BaseCommand):
    help = "Import inventory data from a CSV file exported from the Excel inventory."

    def add_arguments(self, parser):
        parser.add_argument(
            "csv_path",
            type=str,
            help="Path to the CSV file.",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Parse the CSV without saving anything to the database.",
        )
        parser.add_argument(
            "--update-existing-stock",
            action="store_true",
            help="Try to update matching stock entries instead of always creating new ones.",
        )

    @transaction.atomic
    def handle(self, *args, **options):
        csv_path = Path(options["csv_path"])
        dry_run = options["dry_run"]
        update_existing_stock = options["update_existing_stock"]

        if not csv_path.exists():
            logger.error(f"CSV file does not exist: {csv_path}")
            raise CommandError(f"CSV file does not exist: {csv_path}")

        logger.info(f"Reading CSV from: {csv_path}")

        call_kwargs = {"csv_path": csv_path}
        accepts_update_existing_stock = (
            "update_existing_stock" in inspect.signature(import_inventory_csv).parameters
        )
        if accepts_update_existing_stock:
            call_kwargs["update_existing_stock"] = update_existing_stock
        elif update_existing_stock:
            logger.warning(
                "Importer does not accept 'update_existing_stock'; option is ignored."
            )

        try:
            summary = import_inventory_csv(**call_kwargs)
        except Exception as exc:
            logger.error(f"Failed to import CSV: {exc}")
            raise CommandError(str(exc)) from exc

        logger.info(
            f"Finished parsing CSV. "
            f"Imported rows: {summary['imported_rows']}. "
            f"Skipped rows: {summary['skipped_rows']}."
        )

        if dry_run:
            logger.info("Dry run completed successfully. Transaction rolled back.")
            transaction.set_rollback(True)
