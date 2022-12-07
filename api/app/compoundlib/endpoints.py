
from drf_auto_endpoint.endpoints import Endpoint
from drf_auto_endpoint.router import register

from .models import (CompoundLibrary, Compound)
from .serializers import CompoundLibrarySerializer
# from .views import CompoundLibraryViewSet
from core.endpoints import DefaultEndpoint


@register
class CompoundLibraryEndpoint(DefaultEndpoint):
    model = CompoundLibrary
    base_serializer = CompoundLibrarySerializer
    filter_fields = ('name',)


@register
class CompoundEndpoint(DefaultEndpoint):
    model = Compound
