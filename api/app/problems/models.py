from django.db import models
from core.models import Plate, Well
from compoundlib.models import CompoundLibrary


class Problem(models.Model):
    EMPTY_WELL = "empty_well"
    TYPE_CHOICES = [
        (EMPTY_WELL, "Empty Well"),
    ]
    type = models.CharField(max_length=100, choices=TYPE_CHOICES)
    well = models.ForeignKey(Well, on_delete=models.CASCADE, null=True, blank=True)
    well_hr = models.CharField(max_length=10, null=True, blank=True)
    plate = models.ForeignKey(Plate, on_delete=models.CASCADE, null=True, blank=True)
    library = models.ForeignKey(
        CompoundLibrary, on_delete=models.CASCADE, null=True, blank=True
    )
    status = models.CharField(max_length=100, default="open")
    details = models.JSONField(null=True, blank=True)
    show = models.BooleanField(default=True)
