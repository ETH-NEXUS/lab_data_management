from rest_framework import viewsets
from django.db.models import Prefetch
from .models import Well, Plate, Measurement, WellWithdrawal, WellCompound
from .serializers import (PlateSerializer, PlateListSerializer, WellSerializer)


class PlateViewSet(viewsets.ModelViewSet):
    def get_serializer_class(self):
        # if self.action == 'list':
        #     return PlateListSerializer
        # else:
        return PlateSerializer

    def get_queryset(self):
        measurements = Prefetch('measurements', queryset=Measurement.objects.all())
        withdrawals = Prefetch('withdrawals', queryset=WellWithdrawal.objects.select_related('target_well').all())
        donors = Prefetch('donors', queryset=WellWithdrawal.objects.select_related('well').all())
        well_compounds = Prefetch('well_compounds', queryset=WellCompound.objects.select_related('compound').all())
        wells = Prefetch(
            'wells',
            queryset=Well.objects
            .select_related('sample')
            .order_by('position')
            .prefetch_related(well_compounds)
            .prefetch_related(withdrawals)
            .prefetch_related(donors)
            .prefetch_related(measurements)
        )
        return Plate.objects.select_related('dimension', 'experiment', 'library').prefetch_related(wells)

    def filter_queryset(self, queryset):
        return super().filter_queryset(queryset)


class WellViewSet(viewsets.ModelViewSet):
    serializer_class = WellSerializer
    queryset = Well.objects.all()
