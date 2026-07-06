from django.db.models import Prefetch, Q
from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.response import Response

from ..dynamic_models import InventoryStock, Order
from ..serializers.dynamic_models_serializers import OrderDetailSerializer, OrderListSerializer


class OrderViewSet(viewsets.ModelViewSet):
    queryset = (
        Order.objects.all()
        .select_related(
            "material",
            "material__brand",
            "material__manufacturer",
            "material__vendor",
            "material__item_type",
            "material__device_type",
            "order_unit",
            "order_unit__unit",
            "who_ordered",
            "project",
        )
        .prefetch_related(
            "material__attributes",
            Prefetch(
                "created_stock_entries",
                queryset=InventoryStock.objects.select_related(
                    "sector",
                    "sector__room",
                ).prefetch_related(
                    "additional_sectors",
                    "additional_sectors__room",
                ),
            ),
        )
        .order_by("-order_date")
    )

    def get_serializer_class(self):
        if self.action == "list":
            return OrderListSerializer
        return OrderDetailSerializer

    def get_queryset(self):
        queryset = super().get_queryset()

        status_value = self.request.query_params.get("status")
        project_id = self.request.query_params.get("project")
        who_ordered_id = self.request.query_params.get("who_ordered")
        material_id = self.request.query_params.get("material")
        search = self.request.query_params.get("search")

        if status_value:
            queryset = queryset.filter(status=status_value)

        if project_id:
            queryset = queryset.filter(project_id=project_id)

        if who_ordered_id:
            queryset = queryset.filter(who_ordered_id=who_ordered_id)

        if material_id:
            queryset = queryset.filter(material_id=material_id)

        if search:
            queryset = queryset.filter(
                Q(material__product_name__icontains=search)
                | Q(notes__icontains=search)
            )

        return queryset

    @action(detail=False, methods=["get"])
    def pending(self, request):
        queryset = self.get_queryset().filter(
            status__in=[Order.STATUS_ORDERED, Order.STATUS_TENTATIVE]
        )
        serializer = OrderListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)

    @action(detail=False, methods=["get"])
    def arrived(self, request):
        queryset = self.get_queryset().filter(status=Order.STATUS_PRODUCT_ARRIVED)
        serializer = OrderListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)
