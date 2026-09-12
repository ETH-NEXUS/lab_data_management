import logging

from django.http import HttpResponse, JsonResponse
from rest_framework import viewsets
from rest_framework.decorators import action, api_view, permission_classes
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response
from rest_framework.views import APIView
from django.db.models import Prefetch
from .serializers import CompoundLibrarySerializer, CompoundSerializer
from .models import CompoundLibrary, Compound
from core.models import Plate
from core.models import Well
from core.models import Threshold, WellWithdrawal
from django.core import management

logger = logging.getLogger(__name__)


class CompoundLibraryViewSet(viewsets.ModelViewSet):
    serializer_class = CompoundLibrarySerializer
    pagination_class = None

    def get_queryset(self):
        plates = Prefetch("plates", queryset=Plate.objects.all().order_by("barcode"))
        return CompoundLibrary.objects.all().prefetch_related(plates)


class CompoundViewSet(viewsets.ModelViewSet):
    serializer_class = CompoundSerializer
    queryset = Compound.objects.all()

    @action(detail=True, methods=["get"])
    def structure(self, request, pk=None):
        return Response({"src": Compound.objects.get(pk=pk).structure_image})


class RedFlagView(APIView):
    """
    Lists the wells that are running low, grouped by library and plate.
    Returned data example:
    {"Library A": {"PLATE-001": ["A01", "B02"]}}
    """

    permission_classes = [IsAuthenticated]

    def get(self, request, *args, **kwargs):
        plates_with_empty_wells_status = Plate.objects.filter(
            status="empty_wells", library__isnull=False
        ).prefetch_related("library")
        res = {}
        for plate in plates_with_empty_wells_status:
            plate_library_name = plate.library.name
            if plate_library_name not in res:
                res[plate_library_name] = {}
            if plate.barcode not in res[plate_library_name]:
                res[plate_library_name][plate.barcode] = []
            empty_wells = Well.objects.filter(plate=plate, status="empty")
            for well in empty_wells:
                res[plate_library_name][plate.barcode].append(well.hr_position)

        return Response(res)


@api_view(["GET"])
@permission_classes([IsAuthenticated])
def recalculate_status(request):
    """
    Marks the wells that are running low again, for every library plate.
    Returned data example:
    {"status": "ok"}
    """
    try:
        management.call_command("find_problems", "mark_empty_wells")
        return Response({"status": "ok"})
    except Exception:
        # The details belong in the server log, not in the response.
        logger.exception("Recalculating the well statuses failed")
        return Response({"error": "Recalculating the well statuses failed"}, status=500)
