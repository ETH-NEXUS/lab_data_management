from django.db.models import Q
from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.response import Response

from ..serializers.static_models_serializers import (
    MaterialMasterDetailSerializer,
    MaterialMasterListSerializer,
    MaterialUnitSerializer,
)
from ..static_models import MaterialMaster, MaterialUnit


class MaterialMasterViewSet(viewsets.ModelViewSet):
    """
    Main material catalog endpoint.

    Useful for:
    - search
    - detail pages
    - filters by item type / device / manufacturer
    """
    queryset = (
        MaterialMaster.objects.all()
        .select_related(
            "brand",
            "manufacturer",
            "vendor",
            "item_type",
            "device_type",
        )
        .prefetch_related(
            "attributes",
            "units__unit",
        )
        .order_by("product_name")
    )

    def get_serializer_class(self):
        if self.action == "list":
            return MaterialMasterListSerializer
        return MaterialMasterDetailSerializer

    def get_queryset(self):
        queryset = super().get_queryset()

        search = self.request.query_params.get("search")
        manufacturer_id = self.request.query_params.get("manufacturer")
        vendor_id = self.request.query_params.get("vendor")
        item_type_id = self.request.query_params.get("item_type")
        device_type_id = self.request.query_params.get("device_type")
        is_active = self.request.query_params.get("is_active")

        if search:
            queryset = queryset.filter(
                Q(product_name__icontains=search)
                | Q(manufacturer_catalog_number__icontains=search)
                | Q(vendor_catalog_number__icontains=search)
                | Q(description__icontains=search)
                | Q(attributes__name__icontains=search)
            ).distinct()

        if manufacturer_id:
            queryset = queryset.filter(manufacturer_id=manufacturer_id)

        if vendor_id:
            queryset = queryset.filter(vendor_id=vendor_id)

        if item_type_id:
            queryset = queryset.filter(item_type_id=item_type_id)

        if device_type_id:
            queryset = queryset.filter(device_type_id=device_type_id)

        if is_active is not None:
            is_active_bool = is_active.lower() in ("true", "1", "yes")
            queryset = queryset.filter(is_active=is_active_bool)

        return queryset

    @action(detail=False, methods=["get"])
    def active(self, request):
        queryset = self.get_queryset().filter(is_active=True)
        serializer = MaterialMasterListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)

    @action(detail=False, methods=["get"])
    def inactive(self, request):
        queryset = self.get_queryset().filter(is_active=False)
        serializer = MaterialMasterListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)

    @action(detail=False, methods=["get"])
    def by_device(self, request):
        device_name = request.query_params.get("device_name")
        queryset = self.get_queryset()

        if device_name:
            queryset = queryset.filter(device_type__name__iexact=device_name)

        serializer = MaterialMasterListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)


class MaterialUnitViewSet(viewsets.ModelViewSet):
    queryset = (
        MaterialUnit.objects.all()
        .select_related("material", "unit")
        .order_by("material__product_name", "unit__name")
    )
    serializer_class = MaterialUnitSerializer

    def get_queryset(self):
        queryset = super().get_queryset()

        material_id = self.request.query_params.get("material")
        is_base_unit = self.request.query_params.get("is_base_unit")
        is_stock_unit = self.request.query_params.get("is_stock_unit")
        is_order_unit = self.request.query_params.get("is_order_unit")

        if material_id:
            queryset = queryset.filter(material_id=material_id)

        if is_base_unit is not None:
            queryset = queryset.filter(is_base_unit=is_base_unit.lower() in ("true", "1", "yes"))

        if is_stock_unit is not None:
            queryset = queryset.filter(is_stock_unit=is_stock_unit.lower() in ("true", "1", "yes"))

        if is_order_unit is not None:
            queryset = queryset.filter(is_order_unit=is_order_unit.lower() in ("true", "1", "yes"))

        return queryset
