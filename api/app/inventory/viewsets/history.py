from django.db.models import Q
from rest_framework import viewsets

from ..history_models import InventoryChangeRecord
from ..pagination import InventoryHistoryPagination
from ..serializers.dynamic_models_serializers import (
    InventoryChangeRecordDetailSerializer,
    InventoryChangeRecordListSerializer,
)


class InventoryChangeRecordViewSet(viewsets.ReadOnlyModelViewSet):
    pagination_class = InventoryHistoryPagination
    queryset = (
        InventoryChangeRecord.objects.all()
        .select_related(
            "performed_by",
            "inventory_stock",
            "inventory_stock__material",
            "inventory_stock__material__brand",
            "inventory_stock__material__manufacturer",
            "inventory_stock__material__vendor",
            "inventory_stock__material__item_type",
            "inventory_stock__material__device_type",
            "inventory_stock__sector",
            "inventory_stock__sector__room",
            "inventory_stock__stock_unit",
            "inventory_stock__stock_unit__unit",
            "order",
            "order__material",
            "order__material__brand",
            "order__material__manufacturer",
            "order__material__vendor",
            "order__material__item_type",
            "order__material__device_type",
            "order__order_unit",
            "order__order_unit__unit",
            "order__who_ordered",
            "order__project",
            "material_usage",
            "material_usage__project",
            "material_usage__experiment",
            "material_usage__inventory_stock",
            "material_usage__inventory_stock__material",
            "material_usage__inventory_stock__material__brand",
            "material_usage__inventory_stock__material__manufacturer",
            "material_usage__inventory_stock__material__vendor",
            "material_usage__inventory_stock__material__item_type",
            "material_usage__inventory_stock__material__device_type",
            "material_usage__inventory_stock__sector",
            "material_usage__inventory_stock__sector__room",
            "material_usage__inventory_stock__stock_unit",
            "material_usage__inventory_stock__stock_unit__unit",
            "material_usage__usage_unit",
            "material_usage__usage_unit__unit",
            "project",
            "experiment",
            "quantity_unit",
            "quantity_unit__unit",
        )
        .prefetch_related(
            "inventory_stock__additional_sectors",
            "inventory_stock__additional_sectors__room",
            "inventory_stock__material__attributes",
            "inventory_stock__material__units__unit",
            "order__material__attributes",
            "material_usage__inventory_stock__additional_sectors",
            "material_usage__inventory_stock__additional_sectors__room",
            "material_usage__inventory_stock__material__attributes",
            "material_usage__inventory_stock__material__units__unit",
        )
        .order_by("-performed_at")
    )

    def get_serializer_class(self):
        if self.action == "list":
            return InventoryChangeRecordListSerializer
        return InventoryChangeRecordDetailSerializer

    def get_queryset(self):
        queryset = super().get_queryset()

        performed_action = self.request.query_params.get("performed_action")
        performed_by_id = self.request.query_params.get("performed_by")
        inventory_stock_id = self.request.query_params.get("inventory_stock")
        order_id = self.request.query_params.get("order")
        material_usage_id = self.request.query_params.get("material_usage")
        project_id = self.request.query_params.get("project")
        experiment_id = self.request.query_params.get("experiment")
        material_id = self.request.query_params.get("material")
        activity_group = self.request.query_params.get("activity_group")
        search = self.request.query_params.get("search")

        if performed_action:
            queryset = queryset.filter(performed_action=performed_action)

        if performed_by_id:
            queryset = queryset.filter(performed_by_id=performed_by_id)

        if inventory_stock_id:
            queryset = queryset.filter(inventory_stock_id=inventory_stock_id)

        if order_id:
            queryset = queryset.filter(order_id=order_id)

        if material_usage_id:
            queryset = queryset.filter(material_usage_id=material_usage_id)

        if project_id:
            queryset = queryset.filter(project_id=project_id)

        if experiment_id:
            queryset = queryset.filter(experiment_id=experiment_id)

        if material_id:
            queryset = queryset.filter(
                Q(inventory_stock__material_id=material_id)
                | Q(order__material_id=material_id)
                | Q(material_usage__inventory_stock__material_id=material_id)
            )

        if activity_group == "check_in_out":
            queryset = queryset.filter(
                Q(performed_action=InventoryChangeRecord.ACTION_STOCK_CREATED)
                | Q(
                    performed_action=InventoryChangeRecord.ACTION_STOCK_UPDATED,
                    quantity_delta__gt=0,
                )
                | Q(performed_action=InventoryChangeRecord.ACTION_USAGE_CREATED)
            )

        if search:
            queryset = queryset.filter(
                Q(notes__icontains=search)
                | Q(performed_action__icontains=search)
                | Q(inventory_stock__material__product_name__icontains=search)
                | Q(order__material__product_name__icontains=search)
                | Q(material_usage__inventory_stock__material__product_name__icontains=search)
                | Q(project__name__icontains=search)
            )

        return queryset
