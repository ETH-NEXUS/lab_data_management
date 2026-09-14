"""
Read-only models on top of the materialized views in the database.
"""

from django.contrib.postgres.fields import ArrayField
from django.db import connection, models


class DictField(models.JSONField):
    def from_db_value(self, value, expression, connection):
        if isinstance(value, dict):
            return value
        return super().from_db_value(value, expression, connection)


class MaterializedViewModel(models.Model):
    @classmethod
    def refresh(self, concurrently=False):
        """Refresh the materialized view"""
        with connection.cursor() as cursor:
            cursor.execute(
                f"REFRESH MATERIALIZED VIEW {'CONCURRENTLY' if concurrently else ''} {{0}}".format(
                    self._meta.db_table
                )
            )

    class Meta:
        abstract = True
        managed = False


class WellDetail(MaterializedViewModel):
    id = models.BigIntegerField(primary_key=True)
    plate_id = models.BigIntegerField()
    type = models.CharField(max_length=50)
    status = models.TextField()
    position = models.IntegerField()
    hr_position = models.CharField(max_length=10)
    initial_amount = models.FloatField(blank=True, null=True)
    withdrawal = models.FloatField(blank=True, null=True)
    amount = models.FloatField(blank=True, null=True)
    compounds = ArrayField(models.TextField(blank=True, null=True))
    measurements = DictField()

    class Meta:
        db_table = "core_welldetail"
        managed = False


class PlateDetail(MaterializedViewModel):
    id = models.BigIntegerField(primary_key=True)
    num_wells = models.IntegerField()
    measurement_labels = ArrayField(models.TextField(blank=True, null=True))
    measurement_timestamps = DictField()
    stats = DictField()
    overall_stats = DictField()

    class Meta:
        db_table = "core_platedetail"
        managed = False


class ExperimentDetail(MaterializedViewModel):
    id = models.BigIntegerField(primary_key=True)
    project_id = models.BigIntegerField(blank=True, null=True)
    measurement_labels = ArrayField(models.TextField(blank=True, null=True))
    measurement_timestamps = DictField()
    stats = DictField()
    overall_stats = DictField()

    class Meta:
        db_table = "core_experimentdetail"
        managed = False
