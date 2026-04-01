from rest_framework.routers import DefaultRouter

from .views import (
    BrandViewSet,
    DeviceTypeViewSet,
    InventoryStockViewSet,
    ItemTypeViewSet,
    ManufacturerViewSet,
    MaterialAttributeViewSet,
    MaterialMasterViewSet,
    MaterialUnitViewSet,
    MaterialUsageViewSet,
    OrderViewSet,
    RoomViewSet,
    SectorViewSet,
    UnitOfMeasureViewSet,
    VendorViewSet,
)

router = DefaultRouter()
router.register(r"inventory/manufacturers", ManufacturerViewSet, basename="inventory-manufacturer")
router.register(r"inventory/brands", BrandViewSet, basename="inventory-brand")
router.register(r"inventory/vendors", VendorViewSet, basename="inventory-vendor")
router.register(r"inventory/device-types", DeviceTypeViewSet, basename="inventory-device-type")
router.register(r"inventory/item-types", ItemTypeViewSet, basename="inventory-item-type")
router.register(r"inventory/material-attributes", MaterialAttributeViewSet, basename="inventory-material-attribute")
router.register(r"inventory/units-of-measure", UnitOfMeasureViewSet, basename="inventory-unit-of-measure")

router.register(r"inventory/materials", MaterialMasterViewSet, basename="inventory-material")
router.register(r"inventory/material-units", MaterialUnitViewSet, basename="inventory-material-unit")

router.register(r"inventory/rooms", RoomViewSet, basename="inventory-room")
router.register(r"inventory/sectors", SectorViewSet, basename="inventory-sector")
router.register(r"inventory/stocks", InventoryStockViewSet, basename="inventory-stock")
router.register(r"inventory/orders", OrderViewSet, basename="inventory-order")
router.register(r"inventory/material-usages", MaterialUsageViewSet, basename="inventory-material-usage")

urlpatterns = router.urls