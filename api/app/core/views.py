from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.response import Response
from django.db.models import Prefetch
from .models import Well, Plate, Measurement
from .serializers import (PlateSerializer, WellSerializer)


class PlateViewSet(viewsets.ModelViewSet):
    serializer_class = PlateSerializer
    # filterset_fields = ('barcode', )
    # ordering_fields = ('barcode', )

    def get_queryset(self):
        measurements = Prefetch('measurements', queryset=Measurement.objects.select_related('well').all())
        wells = Prefetch('wells', queryset=Well.objects.select_related('sample', 'compound').order_by(
            'position').prefetch_related('from_well').prefetch_related(measurements))
        return Plate.objects.select_related('dimension', 'project', 'library').prefetch_related(wells)

    def filter_queryset(self, queryset):
        return super().filter_queryset(queryset)


class WellViewSet(viewsets.ModelViewSet):
    serializer_class = WellSerializer
    queryset = Well.objects.all()

    @action(detail=True, methods=['get'])
    def structure(self, request, pk=None):
        return Response({
            'src': Well.objects.get(pk=pk).compound.structure_image
        })
