from django.db import models

from core.basemodels import TimeTrackedModel


class PlateTemplateCategory(models.Model):
    name = models.CharField(max_length=50, unique=True)

    def __str__(self):
        return self.name


class PlateTemplate(TimeTrackedModel):
    related_name = "templates"
    name = models.CharField(max_length=50, db_index=True)
    category = models.ForeignKey(
        PlateTemplateCategory,
        on_delete=models.CASCADE,
        related_name=related_name,
    )

    def __str__(self):
        return self.name

    class Meta:
        unique_together = ("name", "category")
