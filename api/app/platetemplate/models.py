from django.db import models

from core.basemodels import TimeTrackedModel


class PlateTemplateCategory(models.Model):
    name = models.CharField(max_length=50)


class PlateTemplate(TimeTrackedModel):
    name = models.CharField(max_length=50)
    category = models.ForeignKey(PlateTemplateCategory, on_delete=models.CASCADE)
