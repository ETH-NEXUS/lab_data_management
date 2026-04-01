from importlib import import_module
from pathlib import Path

from django.core.management.base import BaseCommand, CommandError

from inventory.static_models import (
    Brand,
    DeviceType,
    ItemType,
    Manufacturer,
    MaterialAttribute,
    UnitOfMeasure,
    Vendor,
)

seed_data_module = import_module("inventory.import.utils.seed_data")
bulk_seed_named_values = seed_data_module.bulk_seed_named_values
load_json_file = seed_data_module.load_json_file


class Command(BaseCommand):
    help = "Seed static inventory lookup values from a JSON file."

    def add_arguments(self, parser):
        parser.add_argument(
            "--file",
            type=str,
            default=None,
            help="Optional path to the JSON seed file.",
        )

    def handle(self, *args, **options):
        default_file = (
            Path(__file__).resolve().parents[2]
            / "import"
            / "data"
            / "inventory_static_seed.json"
        )
        file_path = Path(options["file"]) if options["file"] else default_file

        if not file_path.exists():
            raise CommandError(f"Seed file does not exist: {file_path}")

        data = load_json_file(file_path)

        created_manufacturers = bulk_seed_named_values(
            Manufacturer, data.get("manufacturers", [])
        )
        created_brands = bulk_seed_named_values(
            Brand, data.get("brands", [])
        )
        created_vendors = bulk_seed_named_values(
            Vendor, data.get("vendors", [])
        )
        created_device_types = bulk_seed_named_values(
            DeviceType, data.get("device_types", [])
        )
        created_item_types = bulk_seed_named_values(
            ItemType, data.get("item_types", [])
        )
        created_attributes = bulk_seed_named_values(
            MaterialAttribute, data.get("material_attributes", [])
        )
        created_units = bulk_seed_named_values(
            UnitOfMeasure, data.get("units_of_measure", [])
        )

        self.stdout.write(f"Created {created_manufacturers} manufacturers.")
        self.stdout.write(f"Created {created_brands} brands.")
        self.stdout.write(f"Created {created_vendors} vendors.")
        self.stdout.write(f"Created {created_device_types} device types.")
        self.stdout.write(f"Created {created_item_types} item types.")
        self.stdout.write(f"Created {created_attributes} material attributes.")
        self.stdout.write(f"Created {created_units} units of measure.")
        self.stdout.write(self.style.SUCCESS("Static inventory data seeded successfully."))
