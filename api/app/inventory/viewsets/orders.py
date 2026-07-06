from django.db.models import Prefetch, Q
from django.db import transaction
from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.response import Response

from ..dynamic_models import InventoryStock, Order
from ..history_models import InventoryChangeRecord
from ..history_utils import record_inventory_action
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

    def _build_order_history_note(self, order):
        """
        Builds one small readable order context string.

        Returned data examples:
        - "PCR Plates ordered for Project A"
        - "Tips ordered without project"
        """
        material_name = order.material.product_name if order.material_id else "Unknown material"
        project_name = order.project.name if order.project_id else "without project"
        return f"{material_name} ordered for {project_name}"

    def perform_create(self, serializer):
        with transaction.atomic():
            order = serializer.save()
            performed_by = self.request.user if self.request.user.is_authenticated else None

            record_inventory_action(
                performed_action=InventoryChangeRecord.ACTION_ORDER_CREATED,
                performed_by=performed_by,
                order=order,
                project=order.project,
                quantity_delta=order.amount,
                quantity_unit=order.order_unit,
                notes=f"Created order for {self._build_order_history_note(order)}.",
            )

    def perform_update(self, serializer):
        previous_order = self.get_object()
        previous_amount = previous_order.amount
        previous_note = self._build_order_history_note(previous_order)

        with transaction.atomic():
            order = serializer.save()
            performed_by = self.request.user if self.request.user.is_authenticated else None

            record_inventory_action(
                performed_action=InventoryChangeRecord.ACTION_ORDER_UPDATED,
                performed_by=performed_by,
                order=order,
                project=order.project,
                quantity_delta=order.amount - previous_amount,
                quantity_unit=order.order_unit,
                notes=(
                    f"Updated order from {previous_note} "
                    f"to {self._build_order_history_note(order)}."
                ),
            )

    def perform_destroy(self, instance):
        performed_by = self.request.user if self.request.user.is_authenticated else None
        order_note = self._build_order_history_note(instance)

        with transaction.atomic():
            record_inventory_action(
                performed_action=InventoryChangeRecord.ACTION_ORDER_DELETED,
                performed_by=performed_by,
                order=instance,
                project=instance.project,
                quantity_delta=instance.amount,
                quantity_unit=instance.order_unit,
                notes=f"Deleted order for {order_note}.",
            )
            instance.delete()

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
