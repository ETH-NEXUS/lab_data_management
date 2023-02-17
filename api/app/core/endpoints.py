

from drf_auto_endpoint.endpoints import Endpoint
from drf_auto_endpoint.router import register

from .models import (
    Plate,
    PlateDimension,
    Well,
    Project,
    WellCompound,
    WellWithdrawal,
    Measurement,
    MeasurementFeature,
    PlateMapping,
    Experiment, BarcodeSpecification)
from .serializers import (PlateSerializer, ExperimentSerializer, ProjectSerializer, PlateMappingSerializer, BarcodeSpecificationSerializer)
from .views import (PlateViewSet, WellViewSet, PlateMappingViewSet, ExperimentViewSet)


class DefaultEndpoint(Endpoint):
    """The default Endpoint"""
    include_str = False

    def get_url(self):
        """ The core endpoint defaults to not include the application name in the apis url. """
        if hasattr(self, 'url') and self.url is not None:
            return self.url

        return '{}'.format(self.model_name.replace('_', '-'))


@register
class PlateEndpoint(DefaultEndpoint):
    model = Plate
    base_serializer = PlateSerializer
    base_viewset = PlateViewSet
    filter_fields = ('barcode', 'library', 'experiment')
    ordering_fields = ('barcode',)


@register
class PlateDimensionEndpoint(DefaultEndpoint):
    model = PlateDimension


@register
class WellEndpoint(DefaultEndpoint):
    base_viewset = WellViewSet
    model = Well
    filter_fields = ('plate__barcode',)


@register
class ProjectEndpoint(DefaultEndpoint):
    base_serializer = ProjectSerializer
    model = Project


@register
class ExperimentEndpoint(DefaultEndpoint):
    base_serializer = ExperimentSerializer
    model = Experiment
    base_viewset = ExperimentViewSet

@register
class BarcodeSpecification(DefaultEndpoint):
    base_serializer = BarcodeSpecificationSerializer
    model = BarcodeSpecification



@register
class WellCompoundEndpoint(DefaultEndpoint):
    model = WellCompound
    filter_fields = ('well__plate__barcode',)


@register
class WellWithdrawalEndpoint(DefaultEndpoint):
    model = WellWithdrawal


@register
class MeasurementFeatureEndpoint(DefaultEndpoint):
    model = MeasurementFeature


@register
class MeasurementEndpoint(DefaultEndpoint):
    model = Measurement


@register
class PlateMappingEndpoint(DefaultEndpoint):
    base_serializer = PlateMappingSerializer
    base_viewset = PlateMappingViewSet
    model = PlateMapping


