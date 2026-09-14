"""
Projects, their experiments and the barcode specifications of an experiment.
"""

from django.contrib.postgres.fields import ArrayField
from django.db import models

from ..basemodels import TimeTrackedModel


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
