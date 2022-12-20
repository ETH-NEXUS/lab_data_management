from django.db import models
from django.db.models import F, Sum
from django.core.validators import MinValueValidator
from django.utils import timezone
from django.conf import settings
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


class Experiment(TimeTrackedModel):
    related_name = 'experiments'
    name = models.CharField(max_length=50)
    project = models.ForeignKey(Project, on_delete=models.RESTRICT, related_name=related_name)

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
    experiment = models.ForeignKey(Experiment, null=True, blank=True, on_delete=models.RESTRICT, related_name=related_name)
    library = models.ForeignKey(CompoundLibrary, null=True, blank=True, on_delete=models.RESTRICT, related_name=related_name)

    class Meta:
        ordering = ('-id',)

    def __str__(self):
        return f"{self.barcode}"

    def copy(self, target: 'Plate', amount: float = 0):
        """Copy a plate. Same as map but 1-to-1"""
        self.map(MappingList.one_to_one(self.dimension.num_wells, amount), target)

    def map(self, mappingList: MappingList, target: 'Plate'):
        """
        Maps this plate to another plate using a mapping list.
        The mapping list is just a idx (old position) -> value (new position) mapping
        """
        wells = self.wells.all().order_by('position')
        for mapping in mappingList:
            from_well = wells[mapping.from_pos]
            well = Well.objects.create(
                position=mapping.to_pos,
                plate=target,
            )
            well.source_wells.add(from_well),
            well.save()
            for compound in from_well.compounds.all():
                # The amount is a suggestion derived from the distribution
                # rate in the source and the total amount of the withdrawal
                from_well_compound = WellCompound.objects.get(well=from_well, compound=compound)
                WellCompound.objects.create(
                    well=well,
                    compound=compound,
                    amount=round(mapping.amount * from_well_compound.amount / from_well.amount, settings.FLOAT_PRECISION) if from_well.amount > 0 else 0
                )
                # We add a withdrawal to the source well
                WellWithdrawal.objects.create(
                    well=from_well,
                    amount=mapping.amount
                )

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
    # A well can contain multiple compounds
    compounds = models.ManyToManyField(
        Compound, through='WellCompound')
    # From multiple source wells
    source_wells = models.ManyToManyField(
        'self', through='WellRelation', symmetrical=False)

    def __str__(self):
        return f"{self.plate.barcode}: {self.hr_position}"

    @property
    def hr_position(self) -> str:
        return self.plate.dimension.hr_position(self.position)

    @property
    def amount(self) -> float:
        """
        Summarizes the compound amounts, subtracts the withdrawals and 
        returns the total amount of compound in this well.
        """
        amount = self.well_compounds.all().aggregate(Sum('amount'))['amount__sum'] or 0
        withdrawal = self.well_withdrawals.all().aggregate(Sum('amount'))['amount__sum'] or 0
        return round(amount - withdrawal, settings.FLOAT_PRECISION)

    class Meta:
        unique_together = ('plate', 'position')


class WellCompound(models.Model):
    """
    This is the representation of a compound in a well.
    """
    related_name = 'well_compounds'
    well = models.ForeignKey(Well, on_delete=models.RESTRICT, related_name=related_name)
    compound = models.ForeignKey(Compound, on_delete=models.RESTRICT, related_name=related_name)
    amount = models.FloatField(default=0, validators=[MinValueValidator(0)])


class WellWithdrawal(TimeTrackedModel):
    """
    This is the representation of a withdrawal from a well compound.
    """
    related_name = 'well_withdrawals'
    well = models.ForeignKey(Well, on_delete=models.RESTRICT, related_name=related_name)
    amount = models.FloatField()


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
