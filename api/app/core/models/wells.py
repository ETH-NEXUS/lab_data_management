"""
Wells of a plate: their type and sample, the compounds in them and the withdrawals.
"""

from django.conf import settings
from django.core.exceptions import ValidationError
from django.core.validators import MinValueValidator
from django.db import models
from django.db.models import Sum
from django.utils.translation import gettext_lazy as _

from compoundlib.models import Compound
from ..basemodels import TimeTrackedModel


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
    # "Plate" is given as a string on purpose: plates.py imports the well models
    # for Plate.map, so importing Plate here would be a circular import.
    # Django resolves the string to core.Plate once all models are loaded.
    plate = models.ForeignKey(
        "Plate",
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
    is_invalid = models.BooleanField(default=False, null=True, blank=True)

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
        Summarizes the compound amounts that were put into this well,
        without subtracting the withdrawals.
        """
        amount = self.well_compounds.all().aggregate(Sum("amount"))["amount__sum"] or 0
        return amount

    @property
    def current_info(self):
        withdrawals = self.withdrawals.all()
        if self.plate.library is None or not withdrawals:
            return None

        last_withdrawal = withdrawals.latest("created_at")
        current_amount = last_withdrawal.current_amount
        current_dmso = last_withdrawal.current_dmso
        return {
            "current_amount": current_amount,
            "current_dmso": current_dmso,
        }


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
    current_amount = models.FloatField(null=True, blank=True)
    current_dmso = models.FloatField(null=True, blank=True)

    def __str__(self):
        return f"{self.well.hr_position} ({self.amount})"
