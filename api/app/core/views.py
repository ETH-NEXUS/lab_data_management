from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.response import Response
from .models import Well
from .serializers import WellSerializer


class WellViewSet(viewsets.ModelViewSet):
    serializer_class = WellSerializer
    queryset = Well.objects.all()

    @action(detail=True, methods=['get'])
    def structure(self, request, pk=None):
        return Response({
            'src': Well.objects.get(pk=pk).compound.structure_image
        })
