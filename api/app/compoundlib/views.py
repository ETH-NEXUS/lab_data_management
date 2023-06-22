from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.response import Response
from django.db.models import Prefetch
from .serializers import CompoundLibrarySerializer, CompoundSerializer
from .models import CompoundLibrary, Compound
from core.models import Plate


class CompoundLibraryViewSet(viewsets.ModelViewSet):
    serializer_class = CompoundLibrarySerializer
    pagination_class = None

    def get_queryset(self):
        plates = Prefetch("plates", queryset=Plate.objects.all().order_by("barcode"))
        for i in CompoundLibrary.objects.all():
            print(i.name)
        return CompoundLibrary.objects.all().prefetch_related(plates)


class CompoundViewSet(viewsets.ModelViewSet):
    serializer_class = CompoundSerializer
    queryset = Compound.objects.all()

    @action(detail=True, methods=["get"])
    def structure(self, request, pk=None):
        return Response({"src": Compound.objects.get(pk=pk).structure_image})
