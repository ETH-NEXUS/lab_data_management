import logging

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
from core.thresholds import threshold_reasons
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
    Lists the wells that are running low, grouped by library and plate, with
    the values the instrument reported and the thresholds they are below.
    A value is null when the instrument never reported it. Archived plates are
    left out: the lab archives plates it no longer uses.
    Returned data example:
    {"Library A": {"PLATE-001": [
        {"position": "A01", "current_amount": 0, "current_dmso": 0,
         "reasons": ["volume", "dmso"]}
    ]}}
    """

    permission_classes = [IsAuthenticated]

    def well_entry(self, well, threshold):
        """
        Describes one marked well for the response.
        Returned data example:
        {"position": "I12", "current_amount": 1.31, "current_dmso": 94.5,
         "reasons": ["volume"]}
        """
        withdrawals = list(well.withdrawals.all())
        last_withdrawal = withdrawals[0] if withdrawals else None
        current_amount = last_withdrawal.current_amount if last_withdrawal else None
        current_dmso = last_withdrawal.current_dmso if last_withdrawal else None

        reasons = threshold_reasons(
            current_amount, current_dmso, threshold.amount, threshold.dmso
        )

        # A plate without a dimension cannot name its wells like "A1". Showing the
        # position number keeps the whole page from failing because of one plate.
        position = well.hr_position if well.plate.dimension else str(well.position)

        return {
            "position": position,
            "current_amount": current_amount,
            "current_dmso": current_dmso,
            "reasons": reasons,
        }

    def get(self, request, *args, **kwargs):
        threshold = Threshold.current()
        # An unset value (null in the database) counts as not archived, as in the
        # navigation tree and on the plate page.
        plates_with_empty_wells_status = (
            Plate.objects.filter(status="empty_wells", library__isnull=False)
            .exclude(archived=True)
            .select_related("library")
        )

        res = {}
        for plate in plates_with_empty_wells_status:
            res.setdefault(plate.library.name, {})[plate.barcode] = []

        # The newest withdrawal carries the values the well last reported, the
        # same one the recalculation judges the well by.
        newest_withdrawals_first = WellWithdrawal.objects.order_by("-created_at")
        marked_wells = (
            Well.objects.filter(
                plate__in=plates_with_empty_wells_status, status="empty"
            )
            .select_related("plate__dimension", "plate__library")
            .prefetch_related(
                Prefetch("withdrawals", queryset=newest_withdrawals_first)
            )
            .order_by("position")
        )
        for well in marked_wells:
            library_name = well.plate.library.name
            res[library_name][well.plate.barcode].append(
                self.well_entry(well, threshold)
            )

        return Response(res)


# POST, because this writes to the database: unlike GET, DRF checks the CSRF
# token for it, so a link or another site cannot start a recalculation.
@api_view(["POST"])
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
