"""
Removes the invisible BOM character from compound names.

Excel saves CSV files with a BOM character at the start. Before the importer
removed it, the first compound of an imported plate file got a name like
"\ufeffDMSO", a second compound next to the real "DMSO".

For every compound with a BOM in its name:
- if a compound with the clean name exists, the wells are moved to it and the
  BOM compound is deleted (a well that has both keeps one entry with both amounts);
- otherwise the compound is renamed.

Run it with --dry-run first, then without it:

    python manage.py remove_bom_from_compound_names --dry-run
"""

from django.core.management.base import BaseCommand
from django.db import transaction

from compoundlib.models import Compound
from core.models import ExperimentDetail, PlateDetail, WellCompound, WellDetail

BOM = "\ufeff"


class Command(BaseCommand):
    help = "Remove the invisible BOM character from compound names."

    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Show what would be changed, without saving anything",
        )

    def handle(self, *args, **options):
        compounds = Compound.objects.filter(name__contains=BOM).order_by("id")
        if not compounds.exists():
            self.stdout.write("No compound names with a BOM character.")
            return

        with transaction.atomic():
            for compound in compounds:
                self.fix_compound(compound)

            if options["dry_run"]:
                # Everything above ran for real, but nothing is kept
                transaction.set_rollback(True)
                self.stdout.write("Dry run: nothing was saved.")
                return

        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
        ExperimentDetail.refresh(concurrently=True)
        self.stdout.write("Done: materialized views refreshed.")

    def fix_compound(self, compound: Compound) -> None:
        clean_name = compound.name.replace(BOM, "")
        shown_name = compound.name.replace(BOM, "<BOM>")
        clean_compound = (
            Compound.objects.filter(name=clean_name)
            .exclude(id=compound.id)
            .order_by("id")
            .first()
        )

        if clean_compound is None:
            compound.name = clean_name
            compound.save()
            self.stdout.write(
                f"Renamed: '{shown_name}' (id {compound.id}) to '{clean_name}'"
            )
            return

        # delete() clears the id, so it is kept for the output
        compound_id = compound.id
        well_count = 0
        for well_compound in compound.well_compounds.all():
            well_count += 1
            existing = WellCompound.objects.filter(
                well=well_compound.well, compound=clean_compound
            ).first()
            if existing:
                existing.amount += well_compound.amount
                existing.save()
                well_compound.delete()
            else:
                well_compound.compound = clean_compound
                well_compound.save()
        compound.delete()
        self.stdout.write(
            f"Merged: '{shown_name}' (id {compound_id}, {well_count} wells) "
            f"into '{clean_name}' (id {clean_compound.id})"
        )
