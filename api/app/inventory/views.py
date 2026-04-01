from datetime import timedelta

from django.db.models import Prefetch, Q
from django.utils import timezone
from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.response import Response

from .dynamic_models import InventoryStock, MaterialUsage, Order, Room, Sector
from .serializers.dynamic_models_serializers import (
    InventoryStockDetailSerializer,
    InventoryStockListSerializer,
    MaterialUsageDetailSerializer,
    MaterialUsageListSerializer,
    OrderDetailSerializer,
    OrderListSerializer,
    RoomSerializer,
    SectorSerializer,
)
from .serializers.static_models_serializers import (
    BrandSerializer,
    DeviceTypeSerializer,
    ItemTypeSerializer,
    ManufacturerSerializer,
    MaterialAttributeSerializer,
    MaterialMasterDetailSerializer,
    MaterialMasterListSerializer,
    MaterialUnitSerializer,
    UnitOfMeasureSerializer,
)
from .static_models import (
    Brand,
    DeviceType,
    ItemType,
    Manufacturer,
    MaterialAttribute,
    MaterialMaster,
    MaterialUnit,
    UnitOfMeasure,
    Vendor,
)
from .serializers.static_models_serializers import VendorSerializer


# =========================================================
# Static lookup viewsets
# =========================================================

class ManufacturerViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Manufacturer.objects.all().order_by("name")
    serializer_class = ManufacturerSerializer


class BrandViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Brand.objects.all().order_by("name")
    serializer_class = BrandSerializer


class VendorViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Vendor.objects.all().order_by("name")
    serializer_class = VendorSerializer


class DeviceTypeViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = DeviceType.objects.all().order_by("name")
    serializer_class = DeviceTypeSerializer


class ItemTypeViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = ItemType.objects.all().order_by("name")
    serializer_class = ItemTypeSerializer


class MaterialAttributeViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = MaterialAttribute.objects.all().order_by("name")
    serializer_class = MaterialAttributeSerializer


class UnitOfMeasureViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = UnitOfMeasure.objects.all().order_by("name")
    serializer_class = UnitOfMeasureSerializer


# =========================================================
# Material master / units
# =========================================================

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


# =========================================================
# Location viewsets
# =========================================================

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


# =========================================================
# Inventory stock viewset
# =========================================================

class InventoryStockViewSet(viewsets.ModelViewSet):
    """
    Main stock endpoint for the UI.

    This is likely the most important endpoint for the first inventory UI.
    """
    queryset = (
        InventoryStock.objects.all()
        .select_related(
            "material",
            "material__brand",
            "material__manufacturer",
            "material__vendor",
            "material__item_type",
            "material__device_type",
            "sector",
            "sector__room",
            "stock_unit",
            "stock_unit__unit",
        )
        .prefetch_related(
            "material__attributes",
            "material__units__unit",
        )
        .order_by("material__product_name", "sector__room__name", "sector__name")
    )

    def get_serializer_class(self):
        if self.action == "list":
            return InventoryStockListSerializer
        return InventoryStockDetailSerializer

    def get_queryset(self):
        queryset = super().get_queryset()

        search = self.request.query_params.get("search")
        room_id = self.request.query_params.get("room")
        sector_id = self.request.query_params.get("sector")
        item_type_id = self.request.query_params.get("item_type")
        device_type_id = self.request.query_params.get("device_type")
        manufacturer_id = self.request.query_params.get("manufacturer")
        vendor_id = self.request.query_params.get("vendor")
        is_favorite = self.request.query_params.get("is_favorite")
        inventory_status = self.request.query_params.get("inventory_status")
        lot_number = self.request.query_params.get("lot_number")

        if search:
            queryset = queryset.filter(
                Q(material__product_name__icontains=search)
                | Q(material__manufacturer_catalog_number__icontains=search)
                | Q(material__vendor_catalog_number__icontains=search)
                | Q(lot_number__icontains=search)
                | Q(notes__icontains=search)
                | Q(sector__name__icontains=search)
                | Q(sector__room__name__icontains=search)
            ).distinct()

        if room_id:
            queryset = queryset.filter(sector__room_id=room_id)

        if sector_id:
            queryset = queryset.filter(sector_id=sector_id)

        if item_type_id:
            queryset = queryset.filter(material__item_type_id=item_type_id)

        if device_type_id:
            queryset = queryset.filter(material__device_type_id=device_type_id)

        if manufacturer_id:
            queryset = queryset.filter(material__manufacturer_id=manufacturer_id)

        if vendor_id:
            queryset = queryset.filter(material__vendor_id=vendor_id)

        if is_favorite is not None:
            queryset = queryset.filter(is_favorite=is_favorite.lower() in ("true", "1", "yes"))

        if inventory_status == "low":
            queryset = queryset.filter(quantity__lt=models.F("minimum_quantity"))

        if inventory_status == "in_stock":
            queryset = queryset.filter(quantity__gte=models.F("minimum_quantity"))

        if lot_number:
            queryset = queryset.filter(lot_number__icontains=lot_number)

        return queryset

    @action(detail=False, methods=["get"])
    def favorites(self, request):
        queryset = self.get_queryset().filter(is_favorite=True)
        serializer = InventoryStockListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)

    @action(detail=False, methods=["get"])
    def low_stock(self, request):
        queryset = self.get_queryset().filter(quantity__lt=models.F("minimum_quantity"))
        serializer = InventoryStockListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)

    @action(detail=False, methods=["get"])
    def expiring_soon(self, request):
        today = timezone.localdate()
        soon_date = today + timedelta(days=30)

        queryset = self.get_queryset().filter(
            expiry_date__isnull=False,
            expiry_date__gte=today,
            expiry_date__lte=soon_date,
        )
        serializer = InventoryStockListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)

    @action(detail=False, methods=["get"])
    def expired(self, request):
        today = timezone.localdate()

        queryset = self.get_queryset().filter(
            expiry_date__isnull=False,
            expiry_date__lt=today,
        )
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


# =========================================================
# Orders
# =========================================================

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


# =========================================================
# Material usage
# =========================================================

class MaterialUsageViewSet(viewsets.ModelViewSet):
    queryset = (
        MaterialUsage.objects.all()
        .select_related(
            "project",
            "experiment",
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
            "usage_unit",
            "usage_unit__unit",
        )
        .prefetch_related(
            "inventory_stock__material__attributes",
            "inventory_stock__material__units__unit",
        )
        .order_by("-used_at")
    )

    def get_serializer_class(self):
        if self.action == "list":
            return MaterialUsageListSerializer
        return MaterialUsageDetailSerializer

    def get_queryset(self):
        queryset = super().get_queryset()

        project_id = self.request.query_params.get("project")
        experiment_id = self.request.query_params.get("experiment")
        material_id = self.request.query_params.get("material")
        search = self.request.query_params.get("search")

        if project_id:
            queryset = queryset.filter(project_id=project_id)

        if experiment_id:
            queryset = queryset.filter(experiment_id=experiment_id)

        if material_id:
            queryset = queryset.filter(inventory_stock__material_id=material_id)

        if search:
            queryset = queryset.filter(
                Q(project__name__icontains=search)
                | Q(inventory_stock__material__product_name__icontains=search)
                | Q(notes__icontains=search)
            )

        return queryset

    @action(detail=False, methods=["get"])
    def by_project(self, request):
        project_id = request.query_params.get("project")
        queryset = self.get_queryset()

        if project_id:
            queryset = queryset.filter(project_id=project_id)

        serializer = MaterialUsageListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)