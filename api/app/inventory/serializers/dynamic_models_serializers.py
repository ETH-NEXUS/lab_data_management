from datetime import timedelta

from django.contrib.auth import get_user_model
from django.utils import timezone
from rest_framework import serializers

from core.models import Project, Experiment
from core.serializers import SimpleProjectSerializer, SimpleExperimentSerializer

from inventory.dynamic_models import (
    Room,
    Sector,
    InventoryStock,
    Order,
    MaterialUsage,
)
from inventory.static_models import MaterialMaster, MaterialUnit
from inventory.serializers.shared_serializers import UserSummarySerializer
from inventory.serializers.static_models_serializers import (
    MaterialMasterListSerializer,
    MaterialMasterDetailSerializer,
    MaterialUnitSummarySerializer,
)

User = get_user_model()


def get_current_date_safe():
    """
    Returns today's date safely for both timezone-aware and naive datetime setups.

    Context:
    - In some deployments `USE_TZ=False`, so `timezone.now()` is naive.
    - `timezone.localdate()` raises `ValueError` for naive datetime values.

    Returned data examples:
    - aware setup: local calendar date in configured timezone
    - naive setup: `timezone.now().date()`
    """
    now_value = timezone.now()

    if timezone.is_naive(now_value):
        return now_value.date()

    return timezone.localdate(now_value)


class RoomSerializer(serializers.ModelSerializer):
    label = serializers.SerializerMethodField()

    class Meta:
        model = Room
        fields = ("id", "name", "label")

    def get_label(self, obj):
        return obj.name


class SectorSerializer(serializers.ModelSerializer):
    room = RoomSerializer(read_only=True)
    room_id = serializers.PrimaryKeyRelatedField(
        source="room",
        queryset=Room.objects.all(),
        write_only=True,
        required=False,
    )
    label = serializers.SerializerMethodField()

    class Meta:
        model = Sector
        fields = (
            "id",
            "name",
            "room",
            "room_id",
            "label",
        )

    def get_label(self, obj):
        return str(obj)


class SectorSummarySerializer(serializers.ModelSerializer):
    room = RoomSerializer(read_only=True)
    label = serializers.SerializerMethodField()

    class Meta:
        model = Sector
        fields = ("id", "name", "room", "label")

    def get_label(self, obj):
        return str(obj)


class InventoryStockListSerializer(serializers.ModelSerializer):
    """
    Compact serializer for stock list/table/card views.
    """
    material = MaterialMasterListSerializer(read_only=True)
    sector = SectorSummarySerializer(read_only=True)
    stock_unit = MaterialUnitSummarySerializer(read_only=True)

    room = serializers.SerializerMethodField()
    inventory_status = serializers.ReadOnlyField()
    quantity_in_base_units = serializers.ReadOnlyField()
    minimum_quantity_in_base_units = serializers.ReadOnlyField()
    location_label = serializers.SerializerMethodField()
    stock_label = serializers.SerializerMethodField()
    is_low_stock = serializers.SerializerMethodField()
    is_expired = serializers.SerializerMethodField()
    is_expiring_soon = serializers.SerializerMethodField()

    class Meta:
        model = InventoryStock
        fields = (
            "id",
            "material",
            "sector",
            "room",
            "stock_unit",
            "quantity",
            "minimum_quantity",
            "inventory_status",
            "quantity_in_base_units",
            "minimum_quantity_in_base_units",
            "lot_number",
            "expiry_date",
            "is_favorite",
            "notes",
            "location_label",
            "stock_label",
            "is_low_stock",
            "is_expired",
            "is_expiring_soon",
            "created_at",
            "updated_at",
        )

    def get_room(self, obj):
        if not obj.sector_id:
            return None
        return RoomSerializer(obj.sector.room).data

    def get_location_label(self, obj):
        if not obj.sector_id:
            return None
        return str(obj.sector)

    def get_stock_label(self, obj):
        unit_name = None
        if obj.stock_unit_id and obj.stock_unit.unit_id:
            unit_name = obj.stock_unit.unit.name

        if unit_name:
            return f"{obj.quantity} {unit_name}"
        return str(obj.quantity)

    def get_is_low_stock(self, obj):
        return obj.inventory_status == "low"

    def get_is_expired(self, obj):
        if not obj.expiry_date:
            return False
        return obj.expiry_date < get_current_date_safe()

    def get_is_expiring_soon(self, obj):
        if not obj.expiry_date:
            return False

        today = get_current_date_safe()
        soon_date = today + timedelta(days=30)

        return today <= obj.expiry_date <= soon_date


class InventoryStockDetailSerializer(serializers.ModelSerializer):
    """
    Detailed serializer for one stock entry.
    """
    material = MaterialMasterDetailSerializer(read_only=True)
    material_id = serializers.PrimaryKeyRelatedField(
        source="material",
        queryset=MaterialMaster.objects.all(),
        write_only=True,
        required=False,
    )

    sector = SectorSummarySerializer(read_only=True)
    sector_id = serializers.PrimaryKeyRelatedField(
        source="sector",
        queryset=Sector.objects.all(),
        write_only=True,
        required=False,
    )

    stock_unit = MaterialUnitSummarySerializer(read_only=True)
    stock_unit_id = serializers.PrimaryKeyRelatedField(
        source="stock_unit",
        queryset=MaterialUnit.objects.all(),
        write_only=True,
        required=False,
    )

    inventory_status = serializers.ReadOnlyField()
    quantity_in_base_units = serializers.ReadOnlyField()
    minimum_quantity_in_base_units = serializers.ReadOnlyField()
    location_label = serializers.SerializerMethodField()
    stock_label = serializers.SerializerMethodField()

    class Meta:
        model = InventoryStock
        fields = (
            "id",
            "material",
            "material_id",
            "sector",
            "sector_id",
            "stock_unit",
            "stock_unit_id",
            "quantity",
            "minimum_quantity",
            "inventory_status",
            "quantity_in_base_units",
            "minimum_quantity_in_base_units",
            "lot_number",
            "expiry_date",
            "is_favorite",
            "notes",
            "location_label",
            "stock_label",
            "created_at",
            "updated_at",
        )
        read_only_fields = (
            "inventory_status",
            "quantity_in_base_units",
            "minimum_quantity_in_base_units",
            "created_at",
            "updated_at",
        )

    def get_location_label(self, obj):
        if not obj.sector_id:
            return None
        return str(obj.sector)

    def get_stock_label(self, obj):
        unit_name = None
        if obj.stock_unit_id and obj.stock_unit.unit_id:
            unit_name = obj.stock_unit.unit.name

        if unit_name:
            return f"{obj.quantity} {unit_name}"
        return str(obj.quantity)


class OrderListSerializer(serializers.ModelSerializer):
    material = MaterialMasterListSerializer(read_only=True)
    order_unit = MaterialUnitSummarySerializer(read_only=True)
    who_ordered = UserSummarySerializer(read_only=True)
    project = SimpleProjectSerializer(read_only=True)

    status_label = serializers.SerializerMethodField()

    class Meta:
        model = Order
        fields = (
            "id",
            "material",
            "order_unit",
            "amount",
            "order_date",
            "status",
            "status_label",
            "who_ordered",
            "project",
            "notes",
            "created_at",
            "updated_at",
        )

    def get_status_label(self, obj):
        return obj.get_status_display()


class OrderDetailSerializer(serializers.ModelSerializer):
    material = MaterialMasterListSerializer(read_only=True)
    material_id = serializers.PrimaryKeyRelatedField(
        source="material",
        queryset=MaterialMaster.objects.all(),
        write_only=True,
        required=False,
    )

    order_unit = MaterialUnitSummarySerializer(read_only=True)
    order_unit_id = serializers.PrimaryKeyRelatedField(
        source="order_unit",
        queryset=MaterialUnit.objects.all(),
        write_only=True,
        required=False,
    )

    who_ordered = UserSummarySerializer(read_only=True)
    who_ordered_id = serializers.PrimaryKeyRelatedField(
        source="who_ordered",
        queryset=User.objects.all(),
        write_only=True,
        required=False,
        allow_null=True,
    )

    project = SimpleProjectSerializer(read_only=True)
    project_id = serializers.PrimaryKeyRelatedField(
        source="project",
        queryset=Project.objects.all(),
        write_only=True,
        required=False,
        allow_null=True,
    )

    status_label = serializers.SerializerMethodField()

    class Meta:
        model = Order
        fields = (
            "id",
            "material",
            "material_id",
            "order_unit",
            "order_unit_id",
            "amount",
            "order_date",
            "status",
            "status_label",
            "who_ordered",
            "who_ordered_id",
            "project",
            "project_id",
            "notes",
            "created_at",
            "updated_at",
        )
        read_only_fields = ("created_at", "updated_at")

    def get_status_label(self, obj):
        return obj.get_status_display()


class MaterialUsageListSerializer(serializers.ModelSerializer):
    project = SimpleProjectSerializer(read_only=True)
    experiment = SimpleExperimentSerializer(read_only=True)
    inventory_stock = InventoryStockListSerializer(read_only=True)
    usage_unit = MaterialUnitSummarySerializer(read_only=True)

    material = serializers.SerializerMethodField()
    quantity_used_in_base_units = serializers.ReadOnlyField()

    class Meta:
        model = MaterialUsage
        fields = (
            "id",
            "project",
            "experiment",
            "inventory_stock",
            "material",
            "usage_unit",
            "quantity_used",
            "quantity_used_in_base_units",
            "used_at",
            "notes",
        )

    def get_material(self, obj):
        if not obj.inventory_stock_id:
            return None
        return MaterialMasterListSerializer(obj.inventory_stock.material).data


class MaterialUsageDetailSerializer(serializers.ModelSerializer):
    project = SimpleProjectSerializer(read_only=True)
    project_id = serializers.PrimaryKeyRelatedField(
        source="project",
        queryset=Project.objects.all(),
        write_only=True,
        required=False,
    )

    experiment = SimpleExperimentSerializer(read_only=True)
    experiment_id = serializers.PrimaryKeyRelatedField(
        source="experiment",
        queryset=Experiment.objects.all(),
        write_only=True,
        required=False,
        allow_null=True,
    )

    inventory_stock = InventoryStockListSerializer(read_only=True)
    inventory_stock_id = serializers.PrimaryKeyRelatedField(
        source="inventory_stock",
        queryset=InventoryStock.objects.all(),
        write_only=True,
        required=False,
    )

    usage_unit = MaterialUnitSummarySerializer(read_only=True)
    usage_unit_id = serializers.PrimaryKeyRelatedField(
        source="usage_unit",
        queryset=MaterialUnit.objects.all(),
        write_only=True,
        required=False,
    )

    material = serializers.SerializerMethodField()
    quantity_used_in_base_units = serializers.ReadOnlyField()

    class Meta:
        model = MaterialUsage
        fields = (
            "id",
            "project",
            "project_id",
            "experiment",
            "experiment_id",
            "inventory_stock",
            "inventory_stock_id",
            "material",
            "usage_unit",
            "usage_unit_id",
            "quantity_used",
            "quantity_used_in_base_units",
            "used_at",
            "notes",
        )
        read_only_fields = (
            "used_at",
            "quantity_used_in_base_units",
        )

    def get_material(self, obj):
        if not obj.inventory_stock_id:
            return None
        return MaterialMasterListSerializer(obj.inventory_stock.material).data
