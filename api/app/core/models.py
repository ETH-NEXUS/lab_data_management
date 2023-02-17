import math

from compoundlib.models import CompoundLibrary, Compound
from django.conf import settings
from django.core.exceptions import ObjectDoesNotExist
from django.core.validators import MinValueValidator
from django.db import models, transaction
from django.db.models import F, Sum, CheckConstraint, Q
from django.utils.translation import gettext_lazy as _
from platetemplate.models import PlateTemplate
from django.contrib.postgres.fields import ArrayField

from .basemodels import TimeTrackedModel
from .mapping import MappingList
from .mapping import PositionMapper


class MappingError(Exception):
    pass


class Project(TimeTrackedModel):
    name = models.CharField(max_length=50, unique=True)
    description = models.TextField(blank=True, null=True)

    def __str__(self):
        return self.name



class Experiment(TimeTrackedModel):
    related_name = 'experiments'
    name = models.CharField(max_length=50)
    description = models.TextField(blank=True, null=True)
    project = models.ForeignKey(Project, on_delete=models.RESTRICT, related_name=related_name)

    def __str__(self):
        return self.name



    class Meta:
        unique_together = ('name', 'project')


class BarcodeSpecification(TimeTrackedModel):
    related_name = 'barcode_specifications'
    prefix = models.CharField(max_length=100)
    number_of_plates = models.IntegerField()
    sides = ArrayField(models.CharField(max_length=20))
    experiment = models.ForeignKey(Experiment, on_delete=models.RESTRICT, related_name=related_name)

    def __str__(self):
        return self.prefix

    class Meta:
        ordering = ['id']

    # def generate_barcodes(self):
    #     return [{
    #         'NorthBarcode': f"{self.prefix}_{i + 1}${self.prefix}_{i + 1}" if 'North' in self.sides else '',
    #         'SouthBarcode': f"{self.prefix}_{i + 1}${self.prefix}_{i + 1}" if 'South' in self.sides else '',
    #         'EastBarcode': f"{self.prefix}_{i + 1}${self.prefix}_{i + 1}" if 'East' in self.sides else '',
    #         'WestBarcode': f"{self.prefix}_{i + 1}${self.prefix}_{i + 1}" if 'West' in self.sides else '',
    #     } for i in range(self.number_of_plates)]



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
        PlateDimension, on_delete=models.RESTRICT, default=None, null=True)
    # A plate can only be a library, experiment or template plate
    experiment = models.ForeignKey(Experiment, null=True, blank=True, on_delete=models.CASCADE,
                                   related_name=related_name)
    library = models.ForeignKey(
        CompoundLibrary,
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name=related_name)
    template = models.OneToOneField(PlateTemplate, null=True, blank=True, on_delete=models.CASCADE,
                                    related_name='plate')

    class Meta:
        ordering = ('-id',)
        # A plate can only be a library, experiment or template plate
        constraints = [
            CheckConstraint(
                check=Q(experiment__isnull=True) & Q(library__isnull=True) & Q(template__isnull=True) |
                      Q(experiment__isnull=False) & Q(library__isnull=True) & Q(template__isnull=True) |
                      Q(experiment__isnull=True) & Q(library__isnull=False) & Q(template__isnull=True) |
                      Q(experiment__isnull=True) & Q(library__isnull=True) & Q(template__isnull=False),
                name='check_only_library_or_experiment_or_template'
            )
        ]

    def __str__(self):
        return f"{self.barcode}"

    @property
    def num_wells(self):
        return self.dimension.num_wells

    def copy(self, target: 'Plate', amount: float = 0):
        """Copy a plate. Same as map but 1-to-1"""
        self.map(MappingList.one_to_one(self.dimension.num_wells, amount), target)

    def map(self, mappingList: MappingList, target: 'Plate'):
        """
        Maps this plate to another plate using a mapping list.
        """
        with transaction.atomic():
            for mapping in mappingList:
                try:
                    from_well = Well.objects.get(plate=self, position=mapping.from_pos)
                except ObjectDoesNotExist:
                    from_well = None
                # We only need to map wells that are not empty
                if from_well:
                    if not target.dimension:
                        raise MappingError(_("Target plate has no dimension assigned"))
                    if mapping.to_pos >= target.num_wells:
                        raise MappingError(_("Target plate too small"))
                    well, created = Well.objects.get_or_create(
                        position=mapping.to_pos,
                        plate=target,
                    )
                    if created:
                        well.save()
                    for compound in from_well.compounds.all():
                        # The amount is a suggestion derived from the distribution
                        # rate in the source and the total amount of the withdrawal
                        from_well_compound = WellCompound.objects.get(well=from_well, compound=compound)
                        amount = round(mapping.amount * from_well_compound.amount / from_well.initial_amount,
                                       settings.FLOAT_PRECISION) if from_well.initial_amount > 0 else 0
                        try:
                            well_compound = WellCompound.objects.get(well=well, compound=compound)
                            well_compound.amount = F('amount') + amount
                            well_compound.save()
                        except ObjectDoesNotExist:
                            WellCompound.objects.create(
                                well=well,
                                compound=compound,
                                amount=amount
                            )
                        try:
                            well_withdrawal = WellWithdrawal.objects.get(well=from_well, target_well=well)
                            # We add a withdrawal to the source well
                            well_withdrawal.amount = F('amount') + mapping.amount
                            well_withdrawal.save()
                        except ObjectDoesNotExist:
                            WellWithdrawal.objects.create(
                                well=from_well,
                                target_well=well,
                                amount=mapping.amount
                            )
    # TODO: implement unmap!!!


class Sample(TimeTrackedModel):
    name = models.CharField(max_length=50, verbose_name='sample')

    def __str__(self):
        return self.name


class WellType(models.Model):
    name = models.CharField(max_length=50, db_index=True)
    description = models.TextField()


class Well(TimeTrackedModel):
    related_name = 'wells'
    plate = models.ForeignKey(
        Plate, on_delete=models.CASCADE, related_name=related_name)
    position = models.PositiveIntegerField(db_index=True)
    sample = models.ForeignKey(
        Sample, null=True, blank=True, on_delete=models.RESTRICT, related_name=related_name)
    # A well can contain multiple compounds
    compounds = models.ManyToManyField(
        Compound, through='WellCompound')
    type = models.ForeignKey(WellType, on_delete=models.RESTRICT, default=1, db_index=True)

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
        withdrawal = self.withdrawals.all().aggregate(Sum('amount'))['amount__sum'] or 0
        return round(amount - withdrawal, settings.FLOAT_PRECISION)

    @property
    def initial_amount(self) -> float:
        """
        Summarizes the compound amounts, subtracts the withdrawals and
        returns the total amount of compound in this well.
        """
        amount = self.well_compounds.all().aggregate(Sum('amount'))['amount__sum'] or 0
        return amount

    class Meta:
        unique_together = ('plate', 'position')


class WellCompound(models.Model):
    """
    This is the representation of a compound in a well.
    """
    related_name = 'well_compounds'
    well = models.ForeignKey(Well, on_delete=models.CASCADE, related_name=related_name)
    compound = models.ForeignKey(Compound, on_delete=models.RESTRICT, related_name=related_name)
    amount = models.FloatField(default=0, validators=[MinValueValidator(0)])

    def __str__(self):
        return f"{self.well.hr_position}: {self.compound.name}"

    class Meta:
        unique_together = ('well', 'compound')


class WellWithdrawal(TimeTrackedModel):
    """
    This is the representation of a withdrawal from a well compound.
    """
    related_name = 'withdrawals'
    well = models.ForeignKey(Well, on_delete=models.CASCADE, related_name=related_name)
    target_well = models.ForeignKey(Well, null=True, on_delete=models.SET_NULL, related_name='donors')
    amount = models.FloatField()

    def __str__(self):
        return f"{self.well.hr_position} ({self.amount})"


class MeasurementFeature(models.Model):
    name = models.CharField(max_length=50, null=True, blank=True, verbose_name='measurement')
    abbrev = models.CharField(max_length=4, null=True, blank=True)
    unit = models.CharField(max_length=10, null=True, blank=True)


class MeasurementMetadata(models.Model):
    data = models.JSONField()


class Measurement(TimeTrackedModel):
    related_name = 'measurements'
    well = models.ForeignKey(Well, on_delete=models.CASCADE, related_name=related_name)
    feature = models.ForeignKey(MeasurementFeature, on_delete=models.RESTRICT, related_name=related_name)
    meta = models.ForeignKey(MeasurementMetadata, on_delete=models.RESTRICT, related_name=related_name, null=True,
                             blank=True)
    value = models.FloatField()
    identifier = models.CharField(max_length=20, null=True, blank=True)

    def __str__(self):
        return f"{self.feature.abbrev}: {self.value}{self.feature.unit}"

    class Meta:
        unique_together = ('well', 'feature')


class PlateMapping(TimeTrackedModel):
    source_plate = models.ForeignKey(Plate, on_delete=models.CASCADE, related_name='mapped_to_plates')
    target_plate = models.ForeignKey(Plate, on_delete=models.CASCADE, related_name='mapped_from_plates')
    from_column = models.CharField(max_length=50, null=True)
    to_column = models.CharField(max_length=50, null=True)
    amount_column = models.CharField(max_length=50, null=True)
    delimiter = models.CharField(max_length=1, default=',', null=True)
    quotechar = models.CharField(max_length=1, default='"', null=True)
    mapping_file = models.FileField(null=True)
    amount = models.FloatField(default=None, validators=[MinValueValidator(0)], null=True)

    def save(self, *args, **kwargs):
        if not self.pk:
            # Object is created: We apply the mapping to the source plate.
            if self.mapping_file:
                self.source_plate.map(
                    MappingList.from_csv(
                        self.mapping_file.name,
                        self.from_column,
                        self.to_column,
                        self.amount_column,
                        self.delimiter,
                        self.quotechar),
                    self.target_plate)
            else:
                self.source_plate.copy(self.target_plate, self.amount)
        # TODO: What do we do on an update or delete??
        super().save(*args, **kwargs)
