from django.db import models
from django.db.models import F
from django.core.validators import MinValueValidator
from django.utils import timezone
from .mapping import PositionMapper
from .mapping import MappingList
from compoundlib.models import CompoundLibrary, Compound
import math


class MappingError(Exception):
    pass


class AutoDateTimeField(models.DateTimeField):
    def pre_save(self, model_instance, add):
        return timezone.now()


class TimeTrackedModel(models.Model):
    created_at = models.DateTimeField(default=timezone.now, editable=False)
    modified_at = AutoDateTimeField(editable=False)

    class Meta:
        abstract = True


class Project(TimeTrackedModel):
    name = models.CharField(max_length=50)

    def __str__(self):
        return self.name


class Location(TimeTrackedModel):
    name = models.CharField(max_length=50, verbose_name='location')

    def __str__(self):
        return self.name


class PlateDimension(models.Model):
    name = models.CharField(max_length=50, verbose_name='plate dimension')
    rows = models.PositiveIntegerField(validators=[MinValueValidator(1)])
    cols = models.PositiveIntegerField(validators=[MinValueValidator(1)])

    @property
    def num_wells(self):
        return self.cols * self.rows

    def position(self, position: str) -> int:
        row, col = PositionMapper.map(position)
        return (row - 1) * self.cols + (col - 1)

    def row_col(self, position: int) -> tuple[int, int]:
        row = math.floor(position / self.cols) + 1
        col = position - (row - 1) * self.cols + 1
        return row, col

    def hr_position(self, position: int) -> str:
        row, col = self.row_col(position)
        return PositionMapper.unmap(row, col)

    def __str__(self):
        return f"{self.name} ({self.cols}x{self.rows})"


class Plate(TimeTrackedModel):
    related_name = 'plates'
    barcode = models.CharField(max_length=50, unique=True, db_index=True)
    dimension = models.ForeignKey(
        PlateDimension, on_delete=models.RESTRICT)
    project = models.ForeignKey(Project, null=True, blank=True, on_delete=models.RESTRICT, related_name=related_name)
    library = models.ForeignKey(CompoundLibrary, null=True, blank=True, on_delete=models.RESTRICT, related_name=related_name)

    def __str__(self):
        return f"{self.barcode}"

    def copy(self, target: 'Plate', amount: float = 0.0):
        """Copy a plate. Same as map but 1-to-1"""
        self.map(MappingList.one_to_one(self.dimension.num_wells), target, amount)

    def map(self, mappingList: MappingList, target: 'Plate', amount: float = 0.0):
        '''
        Maps this plate to another plate using a mapping list.
        The mapping list is just a idx (old position) -> value (new position) mapping
        '''
        wells = self.wells.all().order_by('position')
        for mapping in mappingList:
            from_well = wells[mapping.from_pos]
            well = Well.objects.create(
                position=mapping.to_pos,
                plate=target,
                compound=from_well.compound,
                amount=amount
            )
            well.from_well.add(from_well),
            well.save()

    # def suggestDimension(self):
    #     num_wells = self.wells.count()
    #     candidates = PlateDimension.objects.filter(num_wells=num_wells)
    #     if candidates.count() > 0:
    #         return candidates.first()
    #     else:
    #         cols = 0
    #         rows = 0
    #         if num_wells % 2 != 0:
    #             if num_wells < 100:
    #                 if num_wells % 8 == 0:
    #                     cols = 8
    #                     rows = num_wells / 8
    #             elif num_wells < 400:
    #                 if num_wells % 16 == 0:
    #                     cols = 16
    #                     rows = num_wells / 16
    #         return PlateDimension.create(name=f"sug_{cols}x{rows}", cols=cols, rows=rows)


class Sample(TimeTrackedModel):
    name = models.CharField(max_length=50, verbose_name='sample')

    def __str__(self):
        return self.name


class Well(TimeTrackedModel):
    related_name = 'wells'
    plate = models.ForeignKey(
        Plate, on_delete=models.RESTRICT, related_name=related_name)
    position = models.PositiveIntegerField()
    sample = models.ForeignKey(
        Sample, null=True, blank=True, on_delete=models.RESTRICT, related_name=related_name)
    compound = models.ForeignKey(
        Compound, null=True, blank=True, on_delete=models.RESTRICT, related_name=related_name)
    amount = models.FloatField(default=0, validators=[MinValueValidator(0)])
    from_well = models.ManyToManyField(
        'self', through='WellRelation', symmetrical=False)

    def __str__(self):
        return f"{self.plate.barcode}: {self.hr_position}"

    @property
    def hr_position(self):
        return self.plate.dimension.hr_position(self.position)

    @property
    def col(self):
        return

    class Meta:
        unique_together = ('plate', 'position')


class Measurement(TimeTrackedModel):
    related_name = 'measurements'
    name = models.CharField(max_length=50, null=True, blank=True, verbose_name='measurement')
    abbrev = models.CharField(max_length=4, null=True, blank=True)
    unit = models.CharField(max_length=10, null=True, blank=True)
    value = models.FloatField()
    well = models.ForeignKey(Well, on_delete=models.RESTRICT, related_name=related_name)

    def __str__(self):
        return f"{self.abbrev}: {self.value}{self.unit}"


class WellRelation(models.Model):
    from_well = models.ForeignKey(
        Well, on_delete=models.RESTRICT, related_name='to_wells')
    to_well = models.ForeignKey(
        Well, on_delete=models.RESTRICT, related_name='from_wells')
    amount = models.FloatField(default=0, validators=[MinValueValidator(0)])

    def __str__(self):
        return f"{self.from_well} -> {self.to_well}"
