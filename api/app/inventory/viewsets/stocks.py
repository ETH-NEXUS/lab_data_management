from datetime import timedelta

from django.db import models
from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response

from ..dynamic_models import InventoryStock, InventoryStockTablePreference
from ..pagination import InventoryStockPagination
from ..serializers.dynamic_models_serializers import (
    InventoryStockDetailSerializer,
    InventoryStockListSerializer,
    InventoryStockTablePreferenceSerializer,
)
from ..stock_queryset import apply_inventory_stock_list_filters, build_inventory_stock_base_queryset
from .shared import get_current_date_safe


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

    def get_queryset(self):
        queryset = build_inventory_stock_base_queryset().filter(is_archived=False)
        return apply_inventory_stock_list_filters(queryset, self.request.query_params)

    @action(detail=False, methods=["get"])
    def favorites(self, request):
        queryset = self.get_queryset().filter(is_favorite=True)
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
        queryset = build_inventory_stock_base_queryset().filter(is_archived=True)
        queryset = apply_inventory_stock_list_filters(queryset, request.query_params)
        page = self.paginate_queryset(queryset)

        if page is not None:
            serializer = InventoryStockListSerializer(page, many=True, context={"request": request})
            return self.get_paginated_response(serializer.data)

        serializer = InventoryStockListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)

    @action(detail=True, methods=["post"])
    def mark_favorite(self, request, pk=None):
        stock = self.get_object()
        stock.is_favorite = True
        stock.save(update_fields=["is_favorite"])
        serializer = InventoryStockDetailSerializer(stock, context={"request": request})
        return Response(serializer.data)

    @action(detail=True, methods=["post"])
    def unmark_favorite(self, request, pk=None):
        stock = self.get_object()
        stock.is_favorite = False
        stock.save(update_fields=["is_favorite"])
        serializer = InventoryStockDetailSerializer(stock, context={"request": request})
        return Response(serializer.data)

    @action(detail=True, methods=["post"])
    def archive(self, request, pk=None):
        stock = self.get_object()
        stock.is_archived = True
        stock.save(update_fields=["is_archived"])
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
