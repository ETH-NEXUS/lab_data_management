

from drf_auto_endpoint.endpoints import Endpoint
from drf_auto_endpoint.router import register

from .models import (Plate, PlateDimension, Well, Project, Experiment, WellCompound, WellWithdrawal, Measurement, MeasurementFeature)
from .serializers import (WellSerializer, PlateSerializer, ExperimentSerializer, ProjectSerializer)
from .views import (PlateViewSet, WellViewSet)


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
    filter_fields = ('barcode',)
    ordering_fields = ('barcode', )


@register
class PlateDimensionEndpoint(DefaultEndpoint):
    model = PlateDimension


@register
class WellEndpoint(DefaultEndpoint):
    base_viewset = WellViewSet
    model = Well
    filter_fields = ('plate__barcode', )


@register
class ProjectEndpoint(DefaultEndpoint):
    base_serializer = ProjectSerializer
    model = Project


@register
class ExperimentEndpoint(DefaultEndpoint):
    base_serializer = ExperimentSerializer
    model = Experiment


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
