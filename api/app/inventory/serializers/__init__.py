from .shared_serializers import (
    SimpleNameSerializer,
    SimpleObjectLabelSerializer,
    UserSummarySerializer,
)

from .static_models_serializers import (
    ManufacturerSerializer,
    BrandSerializer,
    VendorSerializer,
    DeviceTypeSerializer,
    ItemTypeSerializer,
    MaterialAttributeSerializer,
    UnitOfMeasureSerializer,
    MaterialUnitSerializer,
    MaterialUnitSummarySerializer,
    MaterialMasterListSerializer,
    MaterialMasterDetailSerializer,
)

from .dynamic_models_serializers import (
    RoomSerializer,
    SectorSerializer,
    SectorSummarySerializer,
    InventoryChangeRecordDetailSerializer,
    InventoryChangeRecordListSerializer,
    InventoryStockListSerializer,
    InventoryStockDetailSerializer,
    OrderListSerializer,
    OrderDetailSerializer,
    MaterialUsageListSerializer,
    MaterialUsageDetailSerializer,
)
