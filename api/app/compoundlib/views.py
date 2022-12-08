from rest_framework import viewsets
from django.db.models import Prefetch
from .serializers import CompoundLibrarySerializer
from .models import CompoundLibrary
from core.models import Plate


class CompoundLibraryViewSet(viewsets.ModelViewSet):
    serializer_class = CompoundLibrarySerializer

    def get_queryset(self):
        plates = Prefetch('plates', queryset=Plate.objects.all().order_by('barcode'))
        return CompoundLibrary.objects.all().prefetch_related(plates)
