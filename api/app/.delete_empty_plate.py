# Deletes the library plate with the empty barcode '' (made by old SDF imports).
# DRY_RUN = True rolls everything back at the end.
from django.contrib.admin.utils import NestedObjects
from django.db import transaction
from core.models import Plate, Well, WellCompound, WellWithdrawal

DRY_RUN = True
ALLOWED = {"core.Plate", "core.Well", "core.WellCompound"}


class RollBack(Exception):
    pass


try:
    with transaction.atomic():
        plate = Plate.objects.select_for_update().filter(barcode="").first()
        if plate is None:
            print("no plate with an empty barcode")
            raise RollBack()

        collector = NestedObjects(using="default")
        collector.collect([plate])
        affected = {m._meta.label: len(o) for m, o in collector.model_objs.items()}
        withdrawals = (
            WellWithdrawal.objects.filter(well__plate=plate).count()
            + WellWithdrawal.objects.filter(target_well__plate=plate).count()
        )
        print("plate", plate.id, "library", plate.library_id, "will delete", affected)
        # Django 4.2 lists a SET_NULL query here even when it matches no rows,
        # so the rows are counted, not the queries.
        set_null_rows = sum(
            rows.count() for queries in collector.field_updates.values() for rows in queries
        )
        print("withdrawals", withdrawals, "SET_NULL rows", set_null_rows,
              "protected", len(collector.protected))
        if set(affected) - ALLOWED or withdrawals or set_null_rows or collector.protected:
            raise RuntimeError("stopped: the plate has more data than wells and compounds")

        wells_before = Well.objects.count()
        well_compounds_before = WellCompound.objects.count()
        plate.delete()
        wells_gone = wells_before - Well.objects.count()
        well_compounds_gone = well_compounds_before - WellCompound.objects.count()
        print("deleted wells", wells_gone, "well compounds", well_compounds_gone)
        if wells_gone != affected.get("core.Well", 0) or well_compounds_gone != affected.get("core.WellCompound", 0):
            raise RuntimeError("stopped: other rows were deleted than expected")

        if DRY_RUN:
            print("dry run: rolled back")
            raise RollBack()
        print("committed")
except RollBack:
    pass
