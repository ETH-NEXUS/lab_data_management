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
    InventoryStockTablePreference,
    Order,
    MaterialUsage,
)
from inventory.history_models import InventoryChangeRecord
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


class InventoryStockSourceOrderSummarySerializer(serializers.ModelSerializer):
    """
    Minimal order summary embedded on stock responses.
    """
    status_label = serializers.SerializerMethodField()

    class Meta:
        model = Order
        fields = (
            "id",
            "order_date",
            "status",
            "status_label",
            "created_at",
        )

    def get_status_label(self, obj):
        return obj.get_status_display()


class OrderCreatedStockSummarySerializer(serializers.ModelSerializer):
    """
    Minimal stock summary embedded on order responses.
    """
    location_label = serializers.SerializerMethodField()

    class Meta:
        model = InventoryStock
        fields = (
            "id",
            "lot_number",
            "expiry_date",
            "location_label",
            "created_at",
        )

    def get_location_label(self, obj):
        return build_inventory_stock_location_label(obj)


def get_inventory_stock_sectors(obj):
    """
    Returns the stock locations as one ordered list without duplicates.

    Returned data examples:
    - `[Sector(id=3, room='C75', name='3.1')]`
    - `[Sector(id=3, room='C75', name='3.1'), Sector(id=4, room='C75', name='3.2')]`
    """
    sector_map = {}

    if obj.sector_id:
        sector_map[obj.sector_id] = obj.sector

    for sector in obj.additional_sectors.all():
        if sector.id == obj.sector_id:
            continue
        sector_map[sector.id] = sector

    return sorted(
        sector_map.values(),
        key=lambda sector: (
            sector.room.name.lower(),
            sector.name.lower(),
            sector.id,
        ),
    )


def build_inventory_stock_location_label(obj):
    """
    Builds one readable location label for one or more sectors.

    Returned data examples:
    - `'C75 / 3.1'`
    - `'C75 / 3.1, 3.2'`
    - `'C41 / 1.1; C75 / 3.1'`
    """
    sectors = get_inventory_stock_sectors(obj)

    if len(sectors) == 0:
        return None

    grouped_sector_names = {}

    for sector in sectors:
        room_name = sector.room.name
        if room_name not in grouped_sector_names:
            grouped_sector_names[room_name] = []
        grouped_sector_names[room_name].append(sector.name)

    room_labels = []

    for room_name in sorted(grouped_sector_names.keys(), key=lambda value: value.lower()):
        sector_names = grouped_sector_names[room_name]
        room_labels.append(f"{room_name} / {', '.join(sector_names)}")

    return "; ".join(room_labels)


class InventoryStockListSerializer(serializers.ModelSerializer):
    """
    Compact serializer for stock list/table/card views.
    """
    material = MaterialMasterListSerializer(read_only=True)
    sector = SectorSummarySerializer(read_only=True)
    sectors = serializers.SerializerMethodField()
    stock_unit = MaterialUnitSummarySerializer(read_only=True)

    room = serializers.SerializerMethodField()
    inventory_status = serializers.ReadOnlyField()
    # Explicit DecimalField (not ReadOnlyField) so these computed values serialize
    # as strings, matching every other decimal quantity in this API (e.g. `quantity`).
    # max_digits is higher than `quantity`'s own 12: this value is quantity * base_units_per_unit,
    # a product of two 12-digit decimals, so its integer part can be up to twice as wide.
    quantity_in_base_units = serializers.DecimalField(max_digits=24, decimal_places=6, read_only=True)
    minimum_quantity_in_base_units = serializers.DecimalField(max_digits=24, decimal_places=6, read_only=True)
    location_label = serializers.SerializerMethodField()
    stock_label = serializers.SerializerMethodField()
    is_low_stock = serializers.SerializerMethodField()
    is_expired = serializers.SerializerMethodField()
    is_expiring_soon = serializers.SerializerMethodField()
    source_order = InventoryStockSourceOrderSummarySerializer(read_only=True)

    class Meta:
        model = InventoryStock
        fields = (
            "id",
            "material",
            "sector",
            "sectors",
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
            "is_archived",
            "archived_at",
            "notes",
            "source_order",
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

    def get_sectors(self, obj):
        return SectorSummarySerializer(get_inventory_stock_sectors(obj), many=True).data

    def get_location_label(self, obj):
        return build_inventory_stock_location_label(obj)

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
    sectors = serializers.SerializerMethodField()
    sector_ids = serializers.PrimaryKeyRelatedField(
        many=True,
        queryset=Sector.objects.select_related("room").all(),
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
    # Explicit DecimalField (not ReadOnlyField) so these computed values serialize
    # as strings, matching every other decimal quantity in this API (e.g. `quantity`).
    # max_digits is higher than `quantity`'s own 12: this value is quantity * base_units_per_unit,
    # a product of two 12-digit decimals, so its integer part can be up to twice as wide.
    quantity_in_base_units = serializers.DecimalField(max_digits=24, decimal_places=6, read_only=True)
    minimum_quantity_in_base_units = serializers.DecimalField(max_digits=24, decimal_places=6, read_only=True)
    location_label = serializers.SerializerMethodField()
    stock_label = serializers.SerializerMethodField()
    source_order = InventoryStockSourceOrderSummarySerializer(read_only=True)
    source_order_id = serializers.PrimaryKeyRelatedField(
        source="source_order",
        queryset=Order.objects.all(),
        write_only=True,
        required=False,
        allow_null=True,
    )

    class Meta:
        model = InventoryStock
        fields = (
            "id",
            "material",
            "material_id",
            "sector",
            "sector_id",
            "sectors",
            "sector_ids",
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
            "is_archived",
            "archived_at",
            "notes",
            "source_order",
            "source_order_id",
            "location_label",
            "stock_label",
            "created_at",
            "updated_at",
        )
        read_only_fields = (
            "inventory_status",
            "quantity_in_base_units",
            "minimum_quantity_in_base_units",
            "archived_at",
            "created_at",
            "updated_at",
        )

    def validate(self, attrs):
        sector_list = attrs.get("sector_ids")
        primary_sector = attrs.get("sector")

        if sector_list is not None:
            if len(sector_list) == 0:
                raise serializers.ValidationError({
                    "sector_ids": "Select at least one sector.",
                })

            unique_sector_ids = {sector.id for sector in sector_list}
            if len(unique_sector_ids) != len(sector_list):
                raise serializers.ValidationError({
                    "sector_ids": "Each sector can only be selected once.",
                })

            room_ids = {sector.room_id for sector in sector_list}
            if len(room_ids) > 1:
                raise serializers.ValidationError({
                    "sector_ids": "All selected sectors must belong to the same room.",
                })

            attrs["sector"] = sector_list[0]
            attrs["_validated_sector_list"] = sector_list
            return attrs

        if primary_sector is not None:
            attrs["_validated_sector_list"] = [primary_sector]
            return attrs

        if self.instance is not None:
            return attrs

        if primary_sector is None:
            raise serializers.ValidationError({
                "sector_id": "Sector is required.",
            })

        return attrs

    def create(self, validated_data):
        sector_list = validated_data.pop("_validated_sector_list", [])
        validated_data.pop("sector_ids", None)
        stock = super().create(validated_data)
        stock.additional_sectors.set([
            sector for sector in sector_list
            if sector.id != stock.sector_id
        ])
        return stock

    def update(self, instance, validated_data):
        sector_list = validated_data.pop("_validated_sector_list", None)
        validated_data.pop("sector_ids", None)
        stock = super().update(instance, validated_data)

        if sector_list is not None:
            stock.additional_sectors.set([
                sector for sector in sector_list
                if sector.id != stock.sector_id
            ])

        return stock

    def get_sectors(self, obj):
        return SectorSummarySerializer(get_inventory_stock_sectors(obj), many=True).data

    def get_location_label(self, obj):
        return build_inventory_stock_location_label(obj)

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
    created_stock_entries = OrderCreatedStockSummarySerializer(many=True, read_only=True)

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
            "created_stock_entries",
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
    created_stock_entries = OrderCreatedStockSummarySerializer(many=True, read_only=True)

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
            "created_stock_entries",
            "created_at",
            "updated_at",
        )
        read_only_fields = ("created_at", "updated_at")

    def get_status_label(self, obj):
        return obj.get_status_display()


class InventoryStockTablePreferenceSerializer(serializers.ModelSerializer):
    """
    Stores one authenticated user's saved stock table state.
    """

    class Meta:
        model = InventoryStockTablePreference
        fields = (
            "id",
            "table_key",
            "preset",
            "sorting",
            "column_filters",
            "column_order",
            "column_visibility",
            "created_at",
            "updated_at",
        )
        read_only_fields = (
            "id",
            "created_at",
            "updated_at",
        )


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


class InventoryChangeRecordListSerializer(serializers.ModelSerializer):
    performed_by = UserSummarySerializer(read_only=True)
    inventory_stock = InventoryStockListSerializer(read_only=True)
    order = OrderListSerializer(read_only=True)
    material_usage = MaterialUsageListSerializer(read_only=True)
    project = SimpleProjectSerializer(read_only=True)
    experiment = SimpleExperimentSerializer(read_only=True)
    quantity_unit = MaterialUnitSummarySerializer(read_only=True)

    class Meta:
        model = InventoryChangeRecord
        fields = (
            "id",
            "performed_action",
            "performed_by",
            "performed_at",
            "inventory_stock",
            "order",
            "material_usage",
            "project",
            "experiment",
            "quantity_delta",
            "quantity_unit",
            "notes",
        )


class InventoryChangeRecordDetailSerializer(serializers.ModelSerializer):
    performed_by = UserSummarySerializer(read_only=True)
    inventory_stock = InventoryStockDetailSerializer(read_only=True)
    order = OrderDetailSerializer(read_only=True)
    material_usage = MaterialUsageDetailSerializer(read_only=True)
    project = SimpleProjectSerializer(read_only=True)
    experiment = SimpleExperimentSerializer(read_only=True)
    quantity_unit = MaterialUnitSummarySerializer(read_only=True)

    class Meta:
        model = InventoryChangeRecord
        fields = (
            "id",
            "performed_action",
            "performed_by",
            "performed_at",
            "inventory_stock",
            "order",
            "material_usage",
            "project",
            "experiment",
            "quantity_delta",
            "quantity_unit",
            "notes",
        )
