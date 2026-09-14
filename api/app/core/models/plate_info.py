"""
Extra information about a plate: library plate, replicate, time point, cell type.
"""

from django.db import models

from .plates import Plate
from .projects import Experiment


class PlateInfo(models.Model):
    plate = models.ForeignKey(Plate, on_delete=models.CASCADE)
    experiment = models.ForeignKey(Experiment, on_delete=models.CASCADE, null=True)
    lib_plate_barcode = models.TextField()
    label = models.TextField()
    replicate = models.TextField()
    measurement_time = models.DateTimeField()
    cell_type = models.TextField()
    condition = models.TextField()
