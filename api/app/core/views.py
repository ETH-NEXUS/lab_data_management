import csv
from uuid import uuid4

from django.core.exceptions import ValidationError
from django.db import IntegrityError
from compoundlib.serializers import (SimpleCompoundLibrarySerializer)
from django.db.models import Prefetch, Q
from rest_framework import viewsets, views, status
from rest_framework.decorators import action
from rest_framework.response import Response

from .models import (Well, Plate, Measurement, WellWithdrawal, WellCompound, PlateMapping, Experiment, BarcodeSpecification, PlateDimension)
from .serializers import (PlateSerializer, WellSerializer, PlateMappingSerializer, SimpleExperimentSerializer,
                          ExperimentSerializer)



class ProjectViewSet(viewsets.ModelViewSet):
  serializer_class = WellSerializer
  queryset = Well.objects.all()

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
        edges.update(
            {edge_key: {'source': node_key(donor.well), 'target': node_key(well), 'label': donor.amount}})
        appendDonors(nodes, edges, donor.well)

    def appendWithdrawals(nodes, edges, well):
      for withdrawal in well.withdrawals.all():
        if (withdrawal.target_well):
          nodes.update({node_key(withdrawal.target_well): {'name': node_name(withdrawal.target_well)}})
          edge_key = str(uuid4())
          edges.update({edge_key: {'source': node_key(well), 'target': node_key(
              withdrawal.target_well), 'label': withdrawal.amount}})
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


class ExperimentViewSet(viewsets.ModelViewSet):
  serializer_class = ExperimentSerializer
  queryset = Experiment.objects.all()

  @action(detail=False, methods=['post'])
  def barcodes(self, request):
    """ Returns an array of barcodes prepared for the downloadCSV function"""
    data = request.data
    experiment = Experiment.objects.get(pk=data['experiment_id'])
    barcode_specification = BarcodeSpecification(prefix=data['prefix'],
                                                 number_of_plates=data['number_of_plates'],
                                                 sides=data['sides'], experiment=experiment)
    barcode_specification.save()
    return Response(status=status.HTTP_200_OK)

  @action(detail=False, methods=['post'])
  def bulk_add_plates(self, request):
    try:
      experiment_id = int(request.data['experiment_id'])
      barcode_specification_id = int(request.data['barcode_specification_id'])
      plate_dimension_id = int(request.data['plate_dimension_id'])
      experiment = Experiment.objects.get(pk=experiment_id)
      barcode_specification = BarcodeSpecification.objects.get(pk=barcode_specification_id)
      plate_dimension = PlateDimension.objects.get(pk=plate_dimension_id)
      number_of_plates = barcode_specification.number_of_plates
    except Experiment.DoesNotExist:
      return Response({'error': 'Experiment not found'}, status=status.HTTP_404_NOT_FOUND)
    except BarcodeSpecification.DoesNotExist:
      return Response({'error': 'Barcode specification not found'}, status=status.HTTP_404_NOT_FOUND)
    except PlateDimension.DoesNotExist:
      return Response({'error': 'Plate dimension not found'}, status=status.HTTP_404_NOT_FOUND)

    plates = []
    for i in range(number_of_plates):
      plate = Plate(barcode=barcode_specification.get_barcode_by_number(i + 1), experiment=experiment,
                    dimension=plate_dimension)
      plates.append(plate)

    try:
      Plate.objects.bulk_create(plates)
    except IntegrityError:
      return Response({'error': 'Could not save plates to database. Probably you have already added the plates '
                                'with these barcode prefix to the current experiment'},
                      status=status.HTTP_500_INTERNAL_SERVER_ERROR)
    except ValidationError:
      return Response({'error': 'Invalid plate data'}, status=status.HTTP_400_BAD_REQUEST)

    return Response({'success': 'Plates added successfully'}, status=status.HTTP_200_OK)
