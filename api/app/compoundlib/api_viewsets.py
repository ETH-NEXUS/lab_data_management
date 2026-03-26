from django_filters.rest_framework import DjangoFilterBackend
from rest_framework.filters import SearchFilter

from .views import CompoundLibraryViewSet, CompoundViewSet


class CompoundLibraryEndpointViewSet(CompoundLibraryViewSet):
    filter_backends = (DjangoFilterBackend,)
    filterset_fields = ("name",)


class CompoundEndpointViewSet(CompoundViewSet):
    filter_backends = (DjangoFilterBackend, SearchFilter)
    filterset_fields = ("name",)
    search_fields = ("name",)
