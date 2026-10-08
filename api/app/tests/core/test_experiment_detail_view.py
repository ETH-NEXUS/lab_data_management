"""
The overall stats of an experiment (core_experimentdetail) list the reads of a
label in order: the n-th entry holds the n-th read of that label of every plate.
Another label of the same plate must not shift them.
"""

from datetime import datetime

from django.test import TestCase

from core.models import (
    Experiment,
    ExperimentDetail,
    Measurement,
    Plate,
    PlateDetail,
    PlateDimension,
    Project,
    Well,
    WellType,
)

FIRST_READ = datetime(2026, 10, 7, 10, 0)
MIDDLE_READ = datetime(2026, 10, 7, 11, 0)
SECOND_READ = datetime(2026, 10, 7, 12, 0)


class ExperimentOverallStatsTest(TestCase):
    def setUp(self):
        project = Project.objects.create(name="P1")
        self.experiment = Experiment.objects.create(name="Screen", project=project)
        self.dimension = PlateDimension.objects.create(name="dim_1x2", rows=1, cols=2)
        self.well_type = WellType.objects.create(name="C", description="compound")

    def plate_with_values(self, barcode, values):
        """values: {(label, measured_at): [value of well 0, value of well 1]}"""
        plate = Plate.objects.create(
            barcode=barcode, dimension=self.dimension, experiment=self.experiment
        )
        wells = [
            Well.objects.create(plate=plate, position=position, type=self.well_type)
            for position in range(2)
        ]
        for (label, measured_at), well_values in values.items():
            for well, value in zip(wells, well_values):
                Measurement.objects.create(
                    well=well, label=label, value=value, measured_at=measured_at
                )

    def test_another_label_does_not_shift_the_reads(self):
        # Plate 1: two reads of Lum, and a corrected label at the same times
        self.plate_with_values(
            "SP_1",
            {
                ("Lum", FIRST_READ): [100, 200],
                ("Lum", SECOND_READ): [300, 400],
                ("Lum_bc_N1_median", FIRST_READ): [-50, 50],
                ("Lum_bc_N1_median", SECOND_READ): [-60, 60],
            },
        )
        # Plate 2: one read of Lum only
        self.plate_with_values("SP_2", {("Lum", FIRST_READ): [10, 1000]})

        ExperimentDetail.refresh()

        stats = ExperimentDetail.objects.get(id=self.experiment.id).overall_stats
        # First read: both plates; second read: plate 1 only
        self.assertEqual([10, 300], stats["Lum"]["min"])
        self.assertEqual([1000, 400], stats["Lum"]["max"])
        self.assertEqual([-50, -60], stats["Lum_bc_N1_median"]["min"])
        self.assertEqual([50, 60], stats["Lum_bc_N1_median"]["max"])

    def test_the_timestamps_are_in_the_order_of_the_values(self):
        # Several well types: their rows are merged by json_merge, which used a set
        self.plate_with_values(
            "SP_1",
            {("Lum", FIRST_READ): [100, 200], ("Lum", SECOND_READ): [300, 400]},
        )
        Well.objects.filter(position=1).update(
            type=WellType.objects.create(name="N1", description="reference")
        )

        PlateDetail.refresh()
        ExperimentDetail.refresh()

        plate = Plate.objects.get(barcode="SP_1")
        for details in [
            PlateDetail.objects.get(id=plate.id),
            ExperimentDetail.objects.get(id=self.experiment.id),
        ]:
            timestamps = details.measurement_timestamps["Lum"]
            hours = [datetime.fromisoformat(text).hour for text in timestamps]
            self.assertEqual(sorted(hours), hours)
            self.assertEqual(2, len(hours))

    def test_a_read_of_only_one_well_type_keeps_the_time_points_in_order(self):
        # Well 0 (C) is read at 10:00 and 12:00, well 1 (N1) also at 11:00
        self.plate_with_values(
            "SP_1",
            {("Lum", FIRST_READ): [100, 200], ("Lum", SECOND_READ): [300, 400]},
        )
        well = Well.objects.get(position=1)
        well.type = WellType.objects.create(name="N1", description="reference")
        well.save()
        Measurement.objects.create(
            well=well, label="Lum", value=250, measured_at=MIDDLE_READ
        )

        PlateDetail.refresh()

        details = PlateDetail.objects.get(id=Plate.objects.get(barcode="SP_1").id)
        hours = [
            datetime.fromisoformat(text).hour
            for text in details.measurement_timestamps["Lum"]
        ]
        self.assertEqual([10, 11, 12], hours)
