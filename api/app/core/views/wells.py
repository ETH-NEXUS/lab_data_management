"""
Single wells of a plate.
"""

from uuid import uuid4
from rest_framework import viewsets, status
from rest_framework.decorators import action
from rest_framework.response import Response
from ..models import (
    Well,
    PlateDetail,
    WellDetail,
)
from ..serializers import (
    WellSerializer,
    ExperimentDetail,
)
from ..archived_plates import ensure_plate_can_be_changed


class WellViewSet(viewsets.ModelViewSet):
    serializer_class = WellSerializer
    queryset = Well.objects.all()

    # Wells of archived plates cannot be created, changed or deleted through the API
    # (see core/archived_plates.py).
    def perform_create(self, serializer):
        ensure_plate_can_be_changed(serializer.validated_data.get("plate"))
        super().perform_create(serializer)

    def perform_update(self, serializer):
        # Both the plate the well is on and a plate it would be moved to.
        ensure_plate_can_be_changed(serializer.instance.plate)
        ensure_plate_can_be_changed(serializer.validated_data.get("plate"))
        super().perform_update(serializer)

    def destroy(self, request, *args, **kwargs):
        well = self.get_object()
        ensure_plate_can_be_changed(well.plate)
        well.delete()

        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
        ExperimentDetail.refresh(concurrently=True)

        return Response(status=status.HTTP_200_OK)

    @action(detail=True, methods=["get"])
    def mark_as_invalid(self, request, pk=None):
        well = self.get_object()
        ensure_plate_can_be_changed(well.plate)
        if well.is_invalid:
            well.is_invalid = False
        else:
            well.is_invalid = True
        well.save()
        # PlateDetail.refresh(concurrently=True)
        # WellDetail.refresh(concurrently=True)
        # ExperimentDetail.refresh(concurrently=True)
        return Response(status=status.HTTP_200_OK)

    @action(detail=True, methods=["get"])
    def chain(self, request, pk=None):
        def node_key(well):
            return f"{well.plate.barcode}_{well.hr_position}"

        def node_name(well):
            return f"{well.plate.barcode}: {well.hr_position}"

        def appendDonors(nodes, edges, well):
            for donor in well.donors.all():
                nodes.update({node_key(donor.well): {"name": node_name(donor.well)}})
                edge_key = str(uuid4())
                edges.update(
                    {
                        edge_key: {
                            "source": node_key(donor.well),
                            "target": node_key(well),
                            "label": donor.amount,
                        }
                    }
                )
                appendDonors(nodes, edges, donor.well)

        def appendWithdrawals(nodes, edges, well):
            for withdrawal in well.withdrawals.all():
                if withdrawal.target_well:
                    nodes.update(
                        {
                            node_key(withdrawal.target_well): {
                                "name": node_name(withdrawal.target_well)
                            }
                        }
                    )
                    edge_key = str(uuid4())
                    edges.update(
                        {
                            edge_key: {
                                "source": node_key(well),
                                "target": node_key(withdrawal.target_well),
                                "label": withdrawal.amount,
                            }
                        }
                    )
                    appendWithdrawals(nodes, edges, withdrawal.target_well)

        well = Well.objects.get(pk=pk)
        nodes = {node_key(well): {"name": node_name(well), "root": True}}
        edges = {}
        appendDonors(nodes, edges, well)
        appendWithdrawals(nodes, edges, well)

        return Response({"nodes": nodes, "edges": edges}, status=status.HTTP_200_OK)
