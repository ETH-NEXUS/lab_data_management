from django.http import HttpResponse, JsonResponse
from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.response import Response
from django.db.models import Prefetch
from .serializers import CompoundLibrarySerializer, CompoundSerializer
from .models import CompoundLibrary, Compound
from core.models import Plate
from django.views import View
from core.models import Well
from core.models import Threshold, WellWithdrawal
from django.core import management


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


from django.http import JsonResponse


class RedFlagView(View):
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

        return JsonResponse(res)


def recalculate_status(request):
    try:
        management.call_command("find_problems", "mark_empty_wells")
        return JsonResponse({"status": "ok"})
    except Exception as e:
        print("EXCEPTION")
        print(e)
        return JsonResponse({"error": str(e)}, status=500)
