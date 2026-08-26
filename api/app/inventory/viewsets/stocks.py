from datetime import timedelta

from django.contrib.auth import get_user_model
from django.db import models
from django.db import transaction
from django.db.models import Prefetch
from django.utils import timezone
from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response

from ..dynamic_models import InventoryStock, InventoryStockTablePreference
from ..history_models import InventoryChangeRecord
from ..history_utils import record_inventory_action
from ..pagination import InventoryStockPagination
from ..serializers.dynamic_models_serializers import (
    InventoryStockDetailSerializer,
    InventoryStockListSerializer,
    InventoryStockTablePreferenceSerializer,
)
from ..stock_queryset import apply_inventory_stock_list_filters, build_inventory_stock_base_queryset
from .shared import get_current_date_safe

User = get_user_model()


class InventoryStockViewSet(viewsets.ModelViewSet):
    """
    Main stock endpoint for the UI.

    This is likely the most important endpoint for the first inventory UI.
    """
    pagination_class = InventoryStockPagination
    queryset = build_inventory_stock_base_queryset().order_by("-created_at")

    def get_serializer_class(self):
        if self.action == "list":
            return InventoryStockListSerializer
        return InventoryStockDetailSerializer

    def _build_stock_history_note(self, stock):
        """
        Builds one small readable stock context string.

        Returned data examples:
        - "PCR Plates at C75 / 3.1"
        - "Tips at C41 / 1.1"
        """
        material_name = stock.material.product_name if stock.material_id else "Unknown material"
        sector_name = str(stock.sector) if stock.sector_id else "Unknown location"
        return f"{material_name} at {sector_name}"

    def perform_create(self, serializer):
        with transaction.atomic():
            stock = serializer.save()
            performed_by = self.request.user if self.request.user.is_authenticated else None
            source_order = stock.source_order
            project = source_order.project if source_order else None

            record_inventory_action(
                performed_action=InventoryChangeRecord.ACTION_STOCK_CREATED,
                performed_by=performed_by,
                inventory_stock=stock,
                order=source_order,
                project=project,
                quantity_delta=stock.quantity,
                quantity_unit=stock.stock_unit,
                notes=f"Created stock entry for {self._build_stock_history_note(stock)}.",
            )

    def perform_update(self, serializer):
        previous_stock = self.get_object()
        previous_quantity = previous_stock.quantity
        previous_notes = self._build_stock_history_note(previous_stock)

        with transaction.atomic():
            stock = serializer.save()
            performed_by = self.request.user if self.request.user.is_authenticated else None
            quantity_delta = stock.quantity - previous_quantity

            record_inventory_action(
                performed_action=InventoryChangeRecord.ACTION_STOCK_UPDATED,
                performed_by=performed_by,
                inventory_stock=stock,
                quantity_delta=quantity_delta,
                quantity_unit=stock.stock_unit,
                notes=(
                    f"Updated stock entry from {previous_notes} "
                    f"to {self._build_stock_history_note(stock)}."
                ),
            )

    def perform_destroy(self, instance):
        performed_by = self.request.user if self.request.user.is_authenticated else None
        stock_note = self._build_stock_history_note(instance)

        with transaction.atomic():
            record_inventory_action(
                performed_action=InventoryChangeRecord.ACTION_STOCK_DELETED,
                performed_by=performed_by,
                inventory_stock=instance,
                quantity_delta=instance.quantity,
                quantity_unit=instance.stock_unit,
                notes=f"Deleted stock entry for {stock_note}.",
            )
            instance.delete()

    def _get_filtered_stock_queryset(self, include_archived=False):
        queryset = build_inventory_stock_base_queryset()

        if not include_archived:
            queryset = queryset.filter(is_archived=False)

        queryset = apply_inventory_stock_list_filters(
            queryset,
            self.request.query_params,
            self.request.user,
        )

        if self.request.user.is_authenticated:
            return queryset.prefetch_related(
                Prefetch(
                    "favorite_users",
                    queryset=User.objects.filter(id=self.request.user.id),
                )
            )

        return queryset

    def get_queryset(self):
        if self.action in {"retrieve", "archive", "restore"}:
            return self._get_filtered_stock_queryset(include_archived=True)

        return self._get_filtered_stock_queryset()

    @action(detail=False, methods=["get"])
    def favorites(self, request):
        queryset = self.get_queryset().filter(favorite_users=request.user)
        page = self.paginate_queryset(queryset)

        if page is not None:
            serializer = InventoryStockListSerializer(page, many=True, context={"request": request})
            return self.get_paginated_response(serializer.data)

        serializer = InventoryStockListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)

    @action(detail=False, methods=["get"])
    def low_stock(self, request):
        queryset = self.get_queryset().filter(quantity__lt=models.F("minimum_quantity"))
        page = self.paginate_queryset(queryset)

        if page is not None:
            serializer = InventoryStockListSerializer(page, many=True, context={"request": request})
            return self.get_paginated_response(serializer.data)

        serializer = InventoryStockListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)

    @action(detail=False, methods=["get"])
    def expiring_soon(self, request):
        today = get_current_date_safe()
        soon_date = today + timedelta(days=30)

        queryset = self.get_queryset().filter(
            expiry_date__isnull=False,
            expiry_date__gte=today,
            expiry_date__lte=soon_date,
        )
        page = self.paginate_queryset(queryset)

        if page is not None:
            serializer = InventoryStockListSerializer(page, many=True, context={"request": request})
            return self.get_paginated_response(serializer.data)

        serializer = InventoryStockListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)

    @action(detail=False, methods=["get"])
    def expired(self, request):
        today = get_current_date_safe()

        queryset = self.get_queryset().filter(
            expiry_date__isnull=False,
            expiry_date__lt=today,
        )
        page = self.paginate_queryset(queryset)

        if page is not None:
            serializer = InventoryStockListSerializer(page, many=True, context={"request": request})
            return self.get_paginated_response(serializer.data)

        serializer = InventoryStockListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)

    @action(detail=False, methods=["get"])
    def archived(self, request):
        queryset = self._get_filtered_stock_queryset(include_archived=True).filter(is_archived=True)
        if not request.query_params.get("ordering"):
            queryset = queryset.order_by("-archived_at", "-id")
        page = self.paginate_queryset(queryset)

        if page is not None:
            serializer = InventoryStockListSerializer(page, many=True, context={"request": request})
            return self.get_paginated_response(serializer.data)

        serializer = InventoryStockListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)

    @action(detail=True, methods=["post"])
    def mark_favorite(self, request, pk=None):
        stock = self.get_object()
        with transaction.atomic():
            stock.favorite_users.add(request.user)
            record_inventory_action(
                performed_action=InventoryChangeRecord.ACTION_STOCK_FAVORITED,
                performed_by=request.user if request.user.is_authenticated else None,
                inventory_stock=stock,
                quantity_delta=None,
                quantity_unit=stock.stock_unit,
                notes=f"Marked stock as favorite for {self._build_stock_history_note(stock)}.",
            )
        serializer = InventoryStockDetailSerializer(stock, context={"request": request})
        return Response(serializer.data)

    @action(detail=True, methods=["post"])
    def unmark_favorite(self, request, pk=None):
        stock = self.get_object()
        with transaction.atomic():
            stock.favorite_users.remove(request.user)
            record_inventory_action(
                performed_action=InventoryChangeRecord.ACTION_STOCK_UNFAVORITED,
                performed_by=request.user if request.user.is_authenticated else None,
                inventory_stock=stock,
                quantity_delta=None,
                quantity_unit=stock.stock_unit,
                notes=f"Removed stock favorite flag for {self._build_stock_history_note(stock)}.",
            )
        serializer = InventoryStockDetailSerializer(stock, context={"request": request})
        return Response(serializer.data)

    @action(detail=True, methods=["post"])
    def archive(self, request, pk=None):
        stock = self.get_object()
        with transaction.atomic():
            stock.is_archived = True
            stock.archived_at = timezone.now()
            stock.save(update_fields=["is_archived", "archived_at"])
            record_inventory_action(
                performed_action=InventoryChangeRecord.ACTION_STOCK_ARCHIVED,
                performed_by=request.user if request.user.is_authenticated else None,
                inventory_stock=stock,
                quantity_delta=None,
                quantity_unit=stock.stock_unit,
                notes=f"Archived stock entry for {self._build_stock_history_note(stock)}.",
            )
        serializer = InventoryStockDetailSerializer(stock, context={"request": request})
        return Response(serializer.data)

    @action(detail=True, methods=["post"])
    def restore(self, request, pk=None):
        stock = self.get_object()
        with transaction.atomic():
            stock.is_archived = False
            stock.archived_at = None
            stock.save(update_fields=["is_archived", "archived_at"])
            record_inventory_action(
                performed_action=InventoryChangeRecord.ACTION_STOCK_RESTORED,
                performed_by=request.user if request.user.is_authenticated else None,
                inventory_stock=stock,
                quantity_delta=None,
                quantity_unit=stock.stock_unit,
                notes=f"Restored stock entry for {self._build_stock_history_note(stock)}.",
            )
        serializer = InventoryStockDetailSerializer(stock, context={"request": request})
        return Response(serializer.data)


class InventoryStockTablePreferenceViewSet(viewsets.GenericViewSet):
    """
    Reads and updates the authenticated user's saved inventory stock table state.
    """

    serializer_class = InventoryStockTablePreferenceSerializer
    permission_classes = [IsAuthenticated]

    def get_queryset(self):
        return InventoryStockTablePreference.objects.filter(user=self.request.user)

    def _get_current_preference(self):
        preference, _created = InventoryStockTablePreference.objects.get_or_create(
            user=self.request.user,
            table_key=InventoryStockTablePreference.TABLE_KEY_INVENTORY_STOCK,
        )
        return preference

    @action(detail=False, methods=["get", "put", "patch"])
    def current(self, request):
        preference = self._get_current_preference()

        if request.method == "GET":
            serializer = self.get_serializer(preference)
            return Response(serializer.data)

        serializer = self.get_serializer(preference, data=request.data, partial=request.method == "PATCH")
        serializer.is_valid(raise_exception=True)
        serializer.save(
            user=request.user,
            table_key=InventoryStockTablePreference.TABLE_KEY_INVENTORY_STOCK,
        )
        return Response(serializer.data)
