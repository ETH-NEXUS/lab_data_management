from rest_framework import viewsets, views, status
from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework.parsers import FileUploadParser
from django.db.models import Prefetch, Q
from .models import (Well, Plate, Measurement, WellWithdrawal, WellCompound, PlateMapping)
from .serializers import (PlateSerializer, WellSerializer, PlateMappingSerializer, SimpleExperimentSerializer)
from compoundlib.serializers import (SimpleCompoundLibrarySerializer)
import csv
from uuid import uuid4


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
            .select_related('sample', 'type')
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
        return Response([{
            'label': plate.barcode,
            'value': plate.id,
            'library': SimpleCompoundLibrarySerializer(plate.library).data if plate.library else None,
            'experiment': SimpleExperimentSerializer(plate.experiment).data if plate.experiment else None
        } for plate in Plate.objects.filter(predicate)])

    def filter_queryset(self, queryset):
        return super().filter_queryset(queryset)


class WellViewSet(viewsets.ModelViewSet):
    serializer_class = WellSerializer
    queryset = Well.objects.all()

    @action(detail=True, methods=['get'])
    def chain(self, request, pk=None):
        def node_key(well):
            return f"{well.plate.barcode}_{well.hr_position}"

        def node_name(well):
            return f"{well.plate.barcode}: {well.hr_position}"

        def appendDonors(nodes, edges, well):
            for donor in well.donors.all():
                nodes.update({node_key(donor.well): {'name': node_name(donor.well)}})
                edge_key = str(uuid4())
                edges.update({edge_key: {'source': node_key(donor.well), 'target': node_key(well), 'label': donor.amount}})
                appendDonors(nodes, edges, donor.well)

        def appendWithdrawals(nodes, edges, well):
            for withdrawal in well.withdrawals.all():
                if (withdrawal.target_well):
                    nodes.update({node_key(withdrawal.target_well): {'name': node_name(withdrawal.target_well)}})
                    edge_key = str(uuid4())
                    edges.update({edge_key: {'source': node_key(well), 'target': node_key(withdrawal.target_well), 'label': withdrawal.amount}})
                    appendWithdrawals(nodes, edges, withdrawal.target_well)

        well = Well.objects.get(pk=pk)
        nodes = {
            node_key(well): {'name': node_name(well), 'root': True}
        }
        edges = {}
        appendDonors(nodes, edges, well)
        appendWithdrawals(nodes, edges, well)

        return Response({'nodes': nodes, 'edges': edges}, status=status.HTTP_200_OK)


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

        return Response(data, status=status.HTTP_200_OK)
