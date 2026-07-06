from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.response import Response

from ..dynamic_models import Room, Sector
from ..serializers.dynamic_models_serializers import RoomSerializer, SectorSerializer


class RoomViewSet(viewsets.ModelViewSet):
    queryset = Room.objects.all().order_by("name")
    serializer_class = RoomSerializer

    @action(detail=True, methods=["get"])
    def sectors(self, request, pk=None):
        room = self.get_object()
        queryset = room.sectors.all().order_by("name")
        serializer = SectorSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)


class SectorViewSet(viewsets.ModelViewSet):
    queryset = Sector.objects.select_related("room").all().order_by("room__name", "name")
    serializer_class = SectorSerializer

    def get_queryset(self):
        queryset = super().get_queryset()
        room_id = self.request.query_params.get("room")

        if room_id:
            queryset = queryset.filter(room_id=room_id)

        return queryset
