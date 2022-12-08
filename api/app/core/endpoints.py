

from drf_auto_endpoint.endpoints import Endpoint
from drf_auto_endpoint.router import register

from .models import (Plate, PlateDimension, Well, WellRelation)
from .serializers import (WellSerializer, PlateSerializer)
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


@register
class WellRelationEndpoint(DefaultEndpoint):
    model = WellRelation
