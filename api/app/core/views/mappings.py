"""
Plate to plate mappings and the preview of a mapping file.
"""

import csv
from rest_framework import viewsets, views, status
from rest_framework.response import Response
from ..models import PlateMapping
from ..serializers import PlateMappingSerializer


class PlateMappingViewSet(viewsets.ModelViewSet):
    serializer_class = PlateMappingSerializer
    queryset = PlateMapping.objects.all()


class MappingPreviewView(views.APIView):
    def post(self, request, format=None):
        delimiter = request.GET.get("delimiter") or ","
        quotechar = request.GET.get("quotechar") or '"'
        for _file in request.data:
            with open(_file, "r") as file:
                reader = csv.DictReader(file, delimiter=delimiter, quotechar=quotechar)
                data = [line for line in reader]
            # We expect only one file!
            break

        return Response(data, status=status.HTTP_200_OK)
