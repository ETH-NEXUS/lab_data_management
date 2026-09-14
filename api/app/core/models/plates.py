"""
Plates, their dimensions and locations, and mapping one plate onto another.
"""

import math
from django.conf import settings
from django.core.exceptions import ObjectDoesNotExist
from django.core.validators import MinValueValidator, RegexValidator
from django.db import models, transaction
from django.db.models import F, CheckConstraint, Q
from django.utils.translation import gettext_lazy as _
from rest_framework.exceptions import APIException

from compoundlib.models import CompoundLibrary
from platetemplate.models import PlateTemplate
from ..basemodels import TimeTrackedModel
from ..utils.plates.mapping import MappingList
from ..utils.plates.positions import PositionMapper
from ..utils.wells.threshold_checks import is_below_threshold
from .materialized_views import PlateDetail, WellDetail
from .projects import Experiment, Project
from .thresholds import Threshold
from .wells import Well, WellCompound, WellType, WellWithdrawal


class MappingError(APIException):
    pass


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

    @classmethod
    def by_num_wells(cls, num_wells: int) -> "PlateDimension":
        if num_wells <= 96:
            return cls.objects.get(name="dim_96_8x12")
        elif num_wells <= 384:
            return cls.objects.get(name="dim_384_16x24")
        elif num_wells <= 1536:
            return cls.objects.get(name="dim_1536_32x48")

    def __str__(self):
        return f"{self.name} ({self.cols}x{self.rows})"


class Plate(TimeTrackedModel):
    related_name = "plates"
    barcode = models.CharField(
        max_length=300,
        unique=True,
        db_index=True,
        validators=[
            RegexValidator(r"[^\s]+", _("Plate barcode must not contain strings!"))
        ],
    )
    dimension = models.ForeignKey(
        PlateDimension, on_delete=models.RESTRICT, default=None, null=True
    )
    # A plate can only be a library, experiment, project or template plate
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
    project = models.ForeignKey(
        Project,
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
    is_control_plate = models.BooleanField(default=False, null=True, blank=True)
    archived = models.BooleanField(default=False, null=True, blank=True)
    status = models.TextField(null=True, blank=True)
    use_as_template_to_select = models.BooleanField(
        default=False, null=True, blank=True
    )

    class Meta:
        ordering = ("-id",)
        # A plate can only be a library, experiment, project or template plate
        constraints = [
            CheckConstraint(
                check=Q(experiment__isnull=True)
                & Q(library__isnull=True)
                & Q(template__isnull=True)
                & Q(project__isnull=True)
                | Q(experiment__isnull=False)
                & Q(library__isnull=True)
                & Q(template__isnull=True)
                & Q(project__isnull=True)
                | Q(experiment__isnull=True)
                & Q(library__isnull=False)
                & Q(template__isnull=True)
                & Q(project__isnull=True)
                | Q(experiment__isnull=True)
                & Q(library__isnull=True)
                & Q(template__isnull=False)
                & Q(project__isnull=True)
                | Q(experiment__isnull=True)
                & Q(library__isnull=True)
                & Q(template__isnull=True)
                & Q(project__isnull=False),
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

    def well_at(self, position: int, create_if_not_exist: bool = False) -> "Well":
        if position > self.num_wells:
            raise ValueError(f"Position must be less than {self.num_wells}.")
        try:
            return self.wells.get(position=position)
        except Well.DoesNotExist:
            if create_if_not_exist:
                return Well.objects.create(plate=self, position=position)
            return None

    def copy(self, target: "Plate", amount: float = 0, map_type: bool = False):
        """Copy a plate. Same as map but 1-to-1"""
        self.map(
            MappingList.one_to_one(self.dimension.num_wells, amount, map_type), target
        )

    def apply_template(self, template_plate: "Plate"):
        """
        Applies a template plate to this plate
        """
        print(f"Applying template {template_plate} to {self}")
        if self.num_wells != template_plate.num_wells:
            raise MappingError(
                f"{_('Template plate must have the same amount of wells')}: {self.num_wells} != {template_plate.num_wells}"
            )
        for position in range(template_plate.num_wells):
            template_well = template_plate.well_at(position)
            well = self.well_at(position, create_if_not_exist=True)
            if well and template_well:
                well.type = template_well.type
                well.save()
            if well and not template_well:
                well.type = WellType.by_name("C")
        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
        return self

    def map(self, mappingList: MappingList, target: "Plate"):
        """
        Maps this plate to another plate using a mapping list.
        """
        thresholds = Threshold.current()
        with transaction.atomic():
            for mapping in mappingList:
                from_well = self.well_at(mapping.from_pos)
                # We only need to map wells that are not empty
                # TODO: If no from_well a destination well could probably be generated anyway..?..
                if from_well:
                    from_well_plate = from_well.plate
                    if not target.dimension:
                        raise MappingError(_("Target plate has no dimension assigned"))
                    if mapping.to_pos >= target.num_wells:
                        raise MappingError(_("Target plate too small"))

                    well, created = Well.objects.update_or_create(
                        position=mapping.to_pos,
                        plate=target,
                        defaults={
                            "status": mapping.status
                        },  # should the well status be taken from mapping status?
                    )

                    if created:  #  To create an id
                        well.save()
                    if mapping.map_type:

                        well.type = from_well.type
                        well.save()
                    for compound in from_well.compounds.all():
                        # The amount is a suggestion derived from the distribution
                        # rate in the source and the total amount of the withdrawal
                        # TODO but for some reason the amount of the source well is always zero
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
                        # TODO: This runs once per compound of the source well,
                        # so a well holding more than one compound would add the
                        # transferred volume to the withdrawal more than once.
                        # Our wells currently hold a single compound each, so it
                        # does not happen yet.
                        try:
                            well_withdrawal = WellWithdrawal.objects.get(
                                well=from_well, target_well=well
                            )
                            # We add a withdrawal to the source well
                            well_withdrawal.amount = F("amount") + mapping.amount
                            # A mapping without an instrument reading (a csv file or a
                            # plate copy) must not erase the last reading of the Echo.
                            if mapping.current_amount is not None:
                                well_withdrawal.current_amount = mapping.current_amount
                            if mapping.current_dmso is not None:
                                well_withdrawal.current_dmso = mapping.current_dmso
                            well_withdrawal.save()
                        except ObjectDoesNotExist:
                            WellWithdrawal.objects.create(
                                well=from_well,
                                target_well=well,
                                amount=mapping.amount,
                                current_amount=mapping.current_amount,
                                current_dmso=mapping.current_dmso,
                            )

                    # A source well that is running low has to be marked no
                    # matter whether this was the first transfer out of it or a
                    # repeated one.
                    if from_well_plate.library and is_below_threshold(
                        mapping.current_amount,
                        mapping.current_dmso,
                        thresholds.amount,
                        thresholds.dmso,
                    ):
                        from_well.status = "empty"
                        from_well.save()
                        # A status someone set by hand, like "disposed", stays: a
                        # plate is only marked when it has no status yet.
                        if not from_well_plate.status:
                            from_well_plate.status = "empty_wells"
                            from_well_plate.save()
            return True
