"""
A mapping from a source plate to a target plate (from a csv file or a plate copy).
"""

from django.core.validators import MinValueValidator
from django.db import models

from ..basemodels import TimeTrackedModel
from .plates import Plate


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
