"""
Measurements of wells, their features and the files they were assigned from.
"""

from django.db import models

from ..basemodels import TimeTrackedModel
from .plates import Plate
from .wells import Well


class MeasurementFeature(models.Model):
    abbrev = models.CharField(max_length=20, unique=True)
    name = models.CharField(
        max_length=50, null=True, blank=True, verbose_name="measurement"
    )
    unit = models.CharField(max_length=10, null=True, blank=True)


class MeasurementAssignment(TimeTrackedModel):
    related_name = "assignments"
    status = models.CharField(max_length=50, default="pending")
    plate = models.ForeignKey(
        Plate, on_delete=models.CASCADE, related_name=related_name
    )
    filename = models.TextField()
    measurement_file = models.FileField(null=True)


class Measurement(TimeTrackedModel):
    related_name = "measurements"
    well = models.ForeignKey(Well, on_delete=models.CASCADE, related_name=related_name)
    feature = models.ForeignKey(
        MeasurementFeature,
        on_delete=models.RESTRICT,
        related_name=related_name,
        null=True,
        blank=True,
    )
    value = models.FloatField()
    label = models.CharField(max_length=50, default="none")
    identifier = models.CharField(max_length=20, null=True, blank=True)
    measured_at = models.DateTimeField(null=True, blank=True)
    measurement_assignment = models.ForeignKey(
        MeasurementAssignment,
        on_delete=models.CASCADE,
        related_name=related_name,
        null=True,
    )

    def __str__(self):
        # if self.feature.abbrev and self.feature.unit:
        #     return f"{self.feature.abbrev}: {self.value}{self.feature.unit}"
        # else:
        return f"{self.label}: {self.value}"

    class Meta:
        unique_together = ("well", "label", "measured_at")
