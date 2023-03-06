from django.db import models
from django.utils import timezone


class AutoDateTimeField(models.DateTimeField):
    def pre_save(self, model_instance, add):
        return timezone.now()


class TimeTrackedModel(models.Model):
    created_at = models.DateTimeField(default=timezone.now, editable=False)
    modified_at = AutoDateTimeField(editable=False)

    class Meta:
        abstract = True
