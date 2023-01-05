from rest_framework import viewsets, views
from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework.parsers import FileUploadParser
from django.db.models import Prefetch, Q
from .models import (Well, Plate, Measurement, WellWithdrawal, WellCompound, PlateMapping)
from .serializers import (PlateSerializer, WellSerializer, PlateMappingSerializer)
import csv


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

    @action(detail=False, methods=['get'])
    def barcodes(self, request):
        """ Returns an array of barcodes"""
        library = request.GET.get('library')
        experiment = request.GET.get('experiment')
        predicate = Q()
        if library:
            predicate |= Q(library__isnull=(library.lower() != 'true'))
        if experiment:
            predicate |= Q(experiment__isnull=(experiment.lower() != 'true'))
        return Response([{'label': plate.barcode, 'value': plate.id} for plate in Plate.objects.filter(predicate)])

    def filter_queryset(self, queryset):
        return super().filter_queryset(queryset)


class WellViewSet(viewsets.ModelViewSet):
    serializer_class = WellSerializer
    queryset = Well.objects.all()


class PlateMappingViewSet(viewsets.ModelViewSet):
    serializer_class = PlateMappingSerializer
    queryset = PlateMapping.objects.all()


class MappingPreviewView(views.APIView):
    def post(self, request, format=None):
        delimiter = request.GET.get('delimiter') or ','
        quotechar = request.GET.get('quotechar') or '"'
        for _file in request.data:
            with open(_file, 'r') as file:
                reader = csv.DictReader(file, delimiter=delimiter, quotechar=quotechar)
                data = [line for line in reader]
            # We expect only one file!
            break

        return Response(data, status=200)
