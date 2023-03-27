import random
from typing import List

import math
import numpy as np
from django.conf import settings
from django.contrib.postgres.fields import ArrayField
from django.core.exceptions import ObjectDoesNotExist
from django.core.validators import MinValueValidator, RegexValidator
from django.db import models, transaction
from django.db.models import F, Sum, CheckConstraint, Q
from django.forms import ValidationError
from django.utils.translation import gettext_lazy as _

from compoundlib.models import CompoundLibrary, Compound
from platetemplate.models import PlateTemplate
from .basemodels import TimeTrackedModel
from .mapping import MappingList
from .mapping import PositionMapper


class MappingError(Exception):
    pass


class Project(TimeTrackedModel):
    name = models.CharField(max_length=50, unique=True)
    description = models.TextField(blank=True, null=True)
    harvest_id = models.IntegerField(blank=True, null=True)
    harvest_notes = models.TextField(blank=True, null=True)

    def __str__(self):
        return self.name


class Experiment(TimeTrackedModel):
    related_name = "experiments"
    name = models.CharField(max_length=50)
    description = models.TextField(blank=True, null=True)
    project = models.ForeignKey(
        Project, on_delete=models.RESTRICT, related_name=related_name
    )

    def __str__(self):
        return self.name

    class Meta:
        unique_together = ("name", "project")

    # it is still not clear what to do if the experiment has several barcode specifications for now the function only
    # checks if the given barcode is in the list of all barcodes we should probably rewrite it at the moment when we
    # will know at what moment of the workflow the aunction will be called and what it should do then
    def check_barcode(self, barcode: str) -> bool:
        barcode_specifications = self.barcode_specifications.all()
        for barcode_specification in barcode_specifications:
            if barcode.startswith(barcode_specification.prefix):
                return True


class BarcodeSpecification(TimeTrackedModel):
    related_name = "barcode_specifications"
    prefix = models.CharField(max_length=100)
    number_of_plates = models.IntegerField(null=True, blank=True)
    sides = ArrayField(models.CharField(max_length=20), null=True, blank=True)
    experiment = models.ForeignKey(
        Experiment, on_delete=models.CASCADE, related_name=related_name
    )

    def __str__(self):
        return self.prefix

    class Meta:
        ordering = ["id"]

    def get_barcode_by_number(self, number: int) -> str:
        return f"{self.prefix}_{number}"


class Location(TimeTrackedModel):
    name = models.CharField(max_length=50, verbose_name="location")

    def __str__(self):
        return self.name


class PlateDimension(models.Model):
    name = models.CharField(max_length=50, verbose_name="plate dimension")
    rows = models.PositiveIntegerField(validators=[MinValueValidator(1)])
    cols = models.PositiveIntegerField(validators=[MinValueValidator(1)])

    @property
    def num_wells(self):
        return self.cols * self.rows

    def position(self, position: str) -> int:
        row, col = PositionMapper.map(position)
        return (row - 1) * self.cols + (col - 1)

    def row_col(self, position: int) -> tuple[int, int]:
        row = math.floor(position / self.cols) + 1
        col = position - (row - 1) * self.cols + 1
        return row, col

    def hr_position(self, position: int) -> str:
        row, col = self.row_col(position)
        return PositionMapper.unmap(row, col)

    def __str__(self):
        return f"{self.name} ({self.cols}x{self.rows})"


class Plate(TimeTrackedModel):
    related_name = "plates"
    barcode = models.CharField(
        max_length=50,
        unique=True,
        db_index=True,
        validators=[
            RegexValidator(r"[^\s]+", _("Plate barcode must not contain strings!"))
        ],
    )
    dimension = models.ForeignKey(
        PlateDimension, on_delete=models.RESTRICT, default=None, null=True
    )
    # A plate can only be a library, experiment or template plate
    experiment = models.ForeignKey(
        Experiment,
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name=related_name,
    )
    library = models.ForeignKey(
        CompoundLibrary,
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name=related_name,
    )
    template = models.OneToOneField(
        PlateTemplate,
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name="plate",
    )

    class Meta:
        ordering = ("-id",)
        # A plate can only be a library, experiment or template plate
        constraints = [
            CheckConstraint(
                check=Q(experiment__isnull=True)
                & Q(library__isnull=True)
                & Q(template__isnull=True)
                | Q(experiment__isnull=False)
                & Q(library__isnull=True)
                & Q(template__isnull=True)
                | Q(experiment__isnull=True)
                & Q(library__isnull=False)
                & Q(template__isnull=True)
                | Q(experiment__isnull=True)
                & Q(library__isnull=True)
                & Q(template__isnull=False),
                name="check_only_library_or_experiment_or_template",
            ),
            models.UniqueConstraint(
                fields=["barcode", "experiment"],
                name="unique_barcode_experiment",
            ),
        ]

    def __str__(self):
        return f"{self.barcode}"

    @property
    def num_wells(self):
        return self.dimension.num_wells

    @property
    def measurements(self):
        return self.__measurements()

    @property
    def z_primes(self):
        results = {}
        for measurement in self.__measurements():
            z_prime = self.z_prime(measurement)
            results[measurement] = z_prime if not math.isnan(z_prime) else "N/A"
        return results

    @property
    def min_max(self):
        results = {}
        for measurement in self.__measurements():
            min = self.min(measurement)
            max = self.max(measurement)
            min_all_types = self.min(measurement, type="All")
            max_all_types = self.max(measurement, type="All")
            results[measurement] = {
                "min": min if not math.isnan(min) else "N/A",
                "max": max if not math.isnan(max) else "N/A",
                "min_all_types": min_all_types
                if not math.isnan(min_all_types)
                else "N/A",
                "max_all_types": max_all_types
                if not math.isnan(max_all_types)
                else "N/A",
            }
        return results

    def well_at(self, position: int, create_if_not_exist: bool = False) -> "Well":
        if position > self.num_wells:
            raise ValueError(f"Position must be less than {self.num_wells}.")
        try:
            return self.wells.get(position=position)
        except Well.DoesNotExist:
            if create_if_not_exist:
                return Well.objects.create(plate=self, position=position)
            return None

    def mean(self, abbrev: str, type: str = "C"):
        measurements = [
            w.measurement(abbrev) for w in self.wells.filter(type__name=type)
        ]
        return np.mean(measurements)

    def std(self, abbrev: str, type: str = "C"):
        measurements = [
            w.measurement(abbrev) for w in self.wells.filter(type__name=type)
        ]
        return np.std(measurements)

    def max(self, abbrev: str, type: str = "C"):
        if type == "All":
            measurements = [
                w.measurement(abbrev)
                for w in self.wells.all()
                if w.measurement(abbrev) is not None
            ]
        else:
            measurements = [
                w.measurement(abbrev) for w in self.wells.filter(type__name=type)
            ]
        return np.max(measurements)

    def min(self, abbrev: str, type: str = "C"):
        if type == "All":
            measurements = [
                w.measurement(abbrev)
                for w in self.wells.all()
                if w.measurement(abbrev) is not None
            ]
        else:
            measurements = [
                w.measurement(abbrev) for w in self.wells.filter(type__name=type)
            ]
        return np.min(measurements)

    def __measurements(self) -> list[str]:
        results = []
        well = self.well_at(random.randrange(self.num_wells))
        for measurement in well.measurements.all():
            if measurement.feature.abbrev not in results:
                results.append(measurement.feature.abbrev)
        return results

    @staticmethod
    def __z_factor(measurement1: List[float], measurement2: List[float]):
        mean1 = np.mean(measurement1)
        mean2 = np.mean(measurement2)
        std1 = np.std(measurement1)
        std2 = np.std(measurement2)
        return 1 - (3 * (std1 + std2)) / abs(mean1 - mean2)

    def z_prime(self, abbrev: str) -> float:
        """
        Calculates the z' (prime) factor of the plate given by a barcode
        """
        measurement_positive = [
            w.measurement(abbrev) for w in self.wells.filter(type__name="P")
        ]
        measurement_negative = [
            w.measurement(abbrev) for w in self.wells.filter(type__name="N")
        ]
        return self.__z_factor(measurement_positive, measurement_negative)

    def z_factor(self, abbrev: str) -> float:
        """
        Calculates the z factor of the plate given by a barcode
        """
        measurement_positive = [
            w.measurement(abbrev) for w in self.wells.filter(type__name="P")
        ]
        measurement_samples = [
            w.measurement(abbrev) for w in self.wells.filter(type__name="C")
        ]
        return self.__z_factor(measurement_positive, measurement_samples)

    def z_scores(self, abbrev: str, type: str = "C") -> List[float]:
        """Returns the z scores of all wells in a list where the index is the position"""
        scores = []
        mean = self.mean(abbrev, type)
        std = self.std(abbrev, type)
        for pos in range(self.num_wells):
            well = self.well_at(pos)
            if well:
                value = (well.measurement(abbrev) - mean) / std
            else:
                value = None
            scores.append(value)
        return scores

    def copy(self, target: "Plate", amount: float = 0):
        """Copy a plate. Same as map but 1-to-1"""
        self.map(MappingList.one_to_one(self.dimension.num_wells, amount), target)

    def apply_template(self, template_plate: "Plate"):
        """
        Applies a template plate to this plate
        """
        if self.num_wells != template_plate.num_wells:
            raise MappingError(
                f"{_('Template plate must have the same amount of wells')}: {self.num_wells} != {template_plate.num_wells}"
            )
        for position in range(template_plate.num_wells):
            template_well = template_plate.well_at(position)
            well = self.well_at(position, create_if_not_exist=True)
            # TODO: At the moment we only map the type
            well.type = template_well.type
            well.save()
        return self

    def map(self, mappingList: MappingList, target: "Plate"):
        """
        Maps this plate to another plate using a mapping list.
        """
        with transaction.atomic():
            for mapping in mappingList:
                from_well = self.well_at(mapping.from_pos)
                # We only need to map wells that are not empty
                if from_well:
                    if not target.dimension:
                        raise MappingError(_("Target plate has no dimension assigned"))
                    if mapping.to_pos >= target.num_wells:
                        raise MappingError(_("Target plate too small"))

                    well, created = Well.objects.update_or_create(
                        position=mapping.to_pos,
                        plate=target,
                        defaults={"status": mapping.status},
                    )
                    if created:
                        # To create an id
                        well.save()
                    for compound in from_well.compounds.all():
                        # The amount is a suggestion derived from the distribution
                        # rate in the source and the total amount of the withdrawal
                        from_well_compound = WellCompound.objects.get(
                            well=from_well, compound=compound
                        )
                        amount = (
                            round(
                                mapping.amount
                                * from_well_compound.amount
                                / from_well.initial_amount,
                                settings.FLOAT_PRECISION,
                            )
                            if from_well.initial_amount > 0
                            else 0
                        )
                        try:
                            well_compound = WellCompound.objects.get(
                                well=well, compound=compound
                            )
                            well_compound.amount = F("amount") + amount
                            well_compound.save()
                        except ObjectDoesNotExist:
                            WellCompound.objects.create(
                                well=well, compound=compound, amount=amount
                            )
                        try:
                            well_withdrawal = WellWithdrawal.objects.get(
                                well=from_well, target_well=well
                            )
                            # We add a withdrawal to the source well
                            well_withdrawal.amount = F("amount") + mapping.amount
                            well_withdrawal.save()
                        except ObjectDoesNotExist:
                            WellWithdrawal.objects.create(
                                well=from_well,
                                target_well=well,
                                amount=mapping.amount,
                            )
            return True  # TODO: implement unmap!!!


class Sample(TimeTrackedModel):
    name = models.CharField(max_length=50, verbose_name="sample")

    def __str__(self):
        return self.name


class WellType(models.Model):
    name = models.CharField(max_length=50, db_index=True)
    description = models.TextField()

    @classmethod
    def by_name(cls, name: str):
        return cls.objects.get(name=name.upper())

    def __str__(self):
        return f"{self.name} ({self.description})"


class Well(TimeTrackedModel):
    related_name = "wells"
    plate = models.ForeignKey(
        Plate,
        on_delete=models.CASCADE,
        related_name=related_name,
        db_index=True,
    )
    position = models.PositiveIntegerField(db_index=True)
    sample = models.ForeignKey(
        Sample,
        null=True,
        blank=True,
        on_delete=models.RESTRICT,
        related_name=related_name,
    )
    # A well can contain multiple compounds
    compounds = models.ManyToManyField(Compound, through="WellCompound")
    type = models.ForeignKey(
        WellType, on_delete=models.RESTRICT, default=1, db_index=True
    )
    status = models.TextField(null=True, blank=True)

    class Meta:
        unique_together = ("plate", "position")

    def clean(self):
        if self.position >= self.plate.num_wells:
            raise ValidationError(
                _(f"Position must not be less than {self.plate.num_wells}.")
            )

    def __str__(self):
        return f"{self.plate.barcode}: {self.hr_position}"

    @property
    def hr_position(self) -> str:
        return self.plate.dimension.hr_position(self.position)

    @property
    def amount(self) -> float:
        """
        Summarizes the compound amounts, subtracts the withdrawals and
        returns the total amount of compound in this well.
        """
        amount = self.well_compounds.all().aggregate(Sum("amount"))["amount__sum"] or 0
        withdrawal = self.withdrawals.all().aggregate(Sum("amount"))["amount__sum"] or 0
        return round(amount - withdrawal, settings.FLOAT_PRECISION)

    @property
    def initial_amount(self) -> float:
        """
        Summarizes the compound amounts, subtracts the withdrawals and
        returns the total amount of compound in this well.
        """
        amount = self.well_compounds.all().aggregate(Sum("amount"))["amount__sum"] or 0
        return amount

    def measurement(self, abbrev: str) -> float:
        """Returns the value of a measurements given by its abbrev"""
        for measurement in self.measurements.all():
            if measurement.feature.abbrev == abbrev:
                return measurement.value

    def z_score(self, abbrev: str, type: str = "C"):
        """Returns the z score for this well"""
        return (
            self.measurement(abbrev) - self.plate.mean(abbrev, type)
        ) / self.plate.std(abbrev, type)


class WellCompound(models.Model):
    """
    This is the representation of a compound in a well.
    """

    related_name = "well_compounds"
    well = models.ForeignKey(Well, on_delete=models.CASCADE, related_name=related_name)
    compound = models.ForeignKey(
        Compound, on_delete=models.RESTRICT, related_name=related_name
    )
    amount = models.FloatField(default=0, validators=[MinValueValidator(0)])

    def __str__(self):
        return f"{self.well.hr_position}: {self.compound.name}"

    class Meta:
        unique_together = ("well", "compound")


class WellWithdrawal(TimeTrackedModel):
    """
    This is the representation of a withdrawal from a well compound.
    """

    related_name = "withdrawals"
    well = models.ForeignKey(Well, on_delete=models.CASCADE, related_name=related_name)
    target_well = models.ForeignKey(
        Well, null=True, on_delete=models.SET_NULL, related_name="donors"
    )
    amount = models.FloatField()

    def __str__(self):
        return f"{self.well.hr_position} ({self.amount})"


class MeasurementFeature(models.Model):
    # TODO: is 20 chars too much for an abbrev
    abbrev = models.CharField(max_length=20, unique=True)
    name = models.CharField(
        max_length=50, null=True, blank=True, verbose_name="measurement"
    )
    unit = models.CharField(max_length=10, null=True, blank=True)


class MeasurementAssignment(models.Model):
    related_name = "assignments"
    assignment_timestamp = models.DateTimeField(auto_now_add=True)
    status = models.CharField(max_length=50, default="pending")
    plate = models.ForeignKey(
        Plate, on_delete=models.CASCADE, related_name=related_name
    )
    filename = models.TextField()
    measurement_file = models.FileField(null=True)
    measurement_timestamp = models.DateTimeField(null=True, blank=True)


class Measurement(TimeTrackedModel):
    related_name = "measurements"
    well = models.ForeignKey(Well, on_delete=models.CASCADE, related_name=related_name)
    feature = models.ForeignKey(
        MeasurementFeature,
        on_delete=models.RESTRICT,
        related_name=related_name,
    )
    value = models.FloatField()
    identifier = models.CharField(max_length=20, null=True, blank=True)
    measurement_assignment = models.ForeignKey(
        MeasurementAssignment,
        on_delete=models.CASCADE,
        related_name=related_name,
        null=True,
    )

    def __str__(self):
        if self.feature.abbrev and self.feature.unit:
            return f"{self.feature.abbrev}: {self.value}{self.feature.unit}"
        else:
            return f"{self.feature.abbrev}: {self.value}"

    class Meta:
        unique_together = ("well", "feature")


class PlateMapping(TimeTrackedModel):
    source_plate = models.ForeignKey(
        Plate, on_delete=models.CASCADE, related_name="mapped_to_plates"
    )
    target_plate = models.ForeignKey(
        Plate, on_delete=models.CASCADE, related_name="mapped_from_plates"
    )
    mapping_file = models.FileField(null=True)

    # If we map from csv
    from_column = models.CharField(max_length=50, null=True)
    to_column = models.CharField(max_length=50, null=True)
    amount_column = models.CharField(max_length=50, null=True)
    delimiter = models.CharField(max_length=1, default=",", null=True)
    quotechar = models.CharField(max_length=1, default='"', null=True)

    # If we copy a plate
    amount = models.FloatField(
        default=None, validators=[MinValueValidator(0)], null=True
    )
    evaluation = models.TextField(null=True, blank=True)
