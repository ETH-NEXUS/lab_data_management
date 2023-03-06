from drf_auto_endpoint.router import register

from .models import CompoundLibrary, Compound
from .serializers import CompoundLibrarySerializer, CompoundSerializer
from .views import CompoundLibraryViewSet, CompoundViewSet
from core.endpoints import DefaultEndpoint


@register
class CompoundLibraryEndpoint(DefaultEndpoint):
    model = CompoundLibrary
    base_serializer = CompoundLibrarySerializer
    base_viewset = CompoundLibraryViewSet
    filter_fields = ("name",)


@register
class CompoundEndpoint(DefaultEndpoint):
    model = Compound
    serializer = CompoundSerializer
    base_viewset = CompoundViewSet
    filter_fields = ("name",)
    search_fields = ("name",)
