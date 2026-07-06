from rest_framework import viewsets

from ..serializers.static_models_serializers import (
    BrandSerializer,
    DeviceTypeSerializer,
    ItemTypeSerializer,
    ManufacturerSerializer,
    MaterialAttributeSerializer,
    UnitOfMeasureSerializer,
    VendorSerializer,
)
from ..static_models import (
    Brand,
    DeviceType,
    ItemType,
    Manufacturer,
    MaterialAttribute,
    UnitOfMeasure,
    Vendor,
)


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
