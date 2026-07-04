from rest_framework import serializers

from inventory.static_models import (
    Manufacturer,
    Brand,
    Vendor,
    DeviceType,
    ItemType,
    MaterialAttribute,
    UnitOfMeasure,
    MaterialMaster,
    MaterialUnit,
)


class ManufacturerSerializer(serializers.ModelSerializer):
    label = serializers.SerializerMethodField()

    class Meta:
        model = Manufacturer
        fields = ("id", "name", "label")

    def get_label(self, obj):
        return str(obj)


class BrandSerializer(serializers.ModelSerializer):
    label = serializers.SerializerMethodField()

    class Meta:
        model = Brand
        fields = ("id", "name", "label")

    def get_label(self, obj):
        return str(obj)


class VendorSerializer(serializers.ModelSerializer):
    label = serializers.SerializerMethodField()

    class Meta:
        model = Vendor
        fields = ("id", "name", "label")

    def get_label(self, obj):
        return str(obj)


class DeviceTypeSerializer(serializers.ModelSerializer):
    label = serializers.SerializerMethodField()

    class Meta:
        model = DeviceType
        fields = ("id", "name", "label")

    def get_label(self, obj):
        return str(obj)


class ItemTypeSerializer(serializers.ModelSerializer):
    label = serializers.SerializerMethodField()

    class Meta:
        model = ItemType
        fields = ("id", "name", "label")

    def get_label(self, obj):
        return str(obj)


class MaterialAttributeSerializer(serializers.ModelSerializer):
    label = serializers.SerializerMethodField()

    class Meta:
        model = MaterialAttribute
        fields = ("id", "name", "label")

    def get_label(self, obj):
        return str(obj)


class UnitOfMeasureSerializer(serializers.ModelSerializer):
    label = serializers.SerializerMethodField()

    class Meta:
        model = UnitOfMeasure
        fields = ("id", "name", "label")

    def get_label(self, obj):
        return str(obj)


class MaterialUnitSerializer(serializers.ModelSerializer):
    """
    Full serializer for material unit definitions.
    """
    unit = UnitOfMeasureSerializer(read_only=True)
    unit_id = serializers.PrimaryKeyRelatedField(
        source="unit",
        queryset=UnitOfMeasure.objects.all(),
        write_only=True,
        required=False,
    )
    display_name = serializers.SerializerMethodField()

    class Meta:
        model = MaterialUnit
        fields = (
            "id",
            "material",
            "unit",
            "unit_id",
            "is_base_unit",
            "is_stock_unit",
            "is_order_unit",
            "base_units_per_unit",
            "notes",
            "display_name",
        )
        read_only_fields = ("material",)

    def get_display_name(self, obj):
        return str(obj)


class MaterialUnitSummarySerializer(serializers.ModelSerializer):
    """
    Smaller serializer for embedding material unit info.
    """
    unit = UnitOfMeasureSerializer(read_only=True)
    display_name = serializers.SerializerMethodField()

    class Meta:
        model = MaterialUnit
        fields = (
            "id",
            "unit",
            "is_base_unit",
            "is_stock_unit",
            "is_order_unit",
            "base_units_per_unit",
            "display_name",
        )

    def get_display_name(self, obj):
        return str(obj)


class MaterialMasterListSerializer(serializers.ModelSerializer):
    """
    Compact serializer for list/table/card views.
    """
    brand = BrandSerializer(read_only=True)
    manufacturer = ManufacturerSerializer(read_only=True)
    vendor = VendorSerializer(read_only=True)
    item_type = ItemTypeSerializer(read_only=True)
    device_type = DeviceTypeSerializer(read_only=True)
    attributes = MaterialAttributeSerializer(read_only=True, many=True)

    capacity_display = serializers.SerializerMethodField()
    label = serializers.SerializerMethodField()
    storage_temperature_label = serializers.SerializerMethodField()

    class Meta:
        model = MaterialMaster
        fields = (
            "id",
            "product_name",
            "label",
            "brand",
            "manufacturer",
            "vendor",
            "manufacturer_catalog_number",
            "vendor_catalog_number",
            "item_type",
            "device_type",
            "attributes",
            "capacity_value",
            "capacity_unit",
            "capacity_display",
            "storage_temperature",
            "storage_temperature_label",
            "safety_data_sheet",
            "default_cost",
            "is_active",
        )

    def get_capacity_display(self, obj):
        if obj.capacity_value is None and not obj.capacity_unit:
            return None
        if obj.capacity_value is None:
            return obj.capacity_unit
        if obj.capacity_unit:
            return f"{obj.capacity_value} {obj.capacity_unit}"
        return str(obj.capacity_value)

    def get_label(self, obj):
        return obj.product_name

    def get_storage_temperature_label(self, obj):
        return obj.get_storage_temperature_display() or None


class MaterialMasterDetailSerializer(serializers.ModelSerializer):
    """
    Detailed serializer for material detail pages.
    """
    brand = BrandSerializer(read_only=True)
    brand_id = serializers.PrimaryKeyRelatedField(
        source="brand",
        queryset=Brand.objects.all(),
        write_only=True,
        required=False,
        allow_null=True,
    )

    manufacturer = ManufacturerSerializer(read_only=True)
    manufacturer_id = serializers.PrimaryKeyRelatedField(
        source="manufacturer",
        queryset=Manufacturer.objects.all(),
        write_only=True,
        required=False,
        allow_null=True,
    )

    vendor = VendorSerializer(read_only=True)
    vendor_id = serializers.PrimaryKeyRelatedField(
        source="vendor",
        queryset=Vendor.objects.all(),
        write_only=True,
        required=False,
        allow_null=True,
    )

    item_type = ItemTypeSerializer(read_only=True)
    item_type_id = serializers.PrimaryKeyRelatedField(
        source="item_type",
        queryset=ItemType.objects.all(),
        write_only=True,
        required=False,
        allow_null=True,
    )

    device_type = DeviceTypeSerializer(read_only=True)
    device_type_id = serializers.PrimaryKeyRelatedField(
        source="device_type",
        queryset=DeviceType.objects.all(),
        write_only=True,
        required=False,
        allow_null=True,
    )

    attributes = MaterialAttributeSerializer(read_only=True, many=True)
    attribute_ids = serializers.PrimaryKeyRelatedField(
        source="attributes",
        queryset=MaterialAttribute.objects.all(),
        many=True,
        write_only=True,
        required=False,
    )

    units = MaterialUnitSummarySerializer(read_only=True, many=True)

    capacity_display = serializers.SerializerMethodField()
    label = serializers.SerializerMethodField()
    storage_temperature_label = serializers.SerializerMethodField()

    class Meta:
        model = MaterialMaster
        fields = (
            "id",
            "product_name",
            "label",
            "brand",
            "brand_id",
            "manufacturer",
            "manufacturer_id",
            "vendor",
            "vendor_id",
            "manufacturer_catalog_number",
            "vendor_catalog_number",
            "item_type",
            "item_type_id",
            "device_type",
            "device_type_id",
            "attributes",
            "attribute_ids",
            "capacity_value",
            "capacity_unit",
            "capacity_display",
            "storage_temperature",
            "storage_temperature_label",
            "safety_data_sheet",
            "description",
            "default_cost",
            "lifetime_days",
            "serial_number",
            "order_number",
            "is_active",
            "units",
        )

    def validate(self, attrs):
        item_type = attrs.get("item_type", getattr(self.instance, "item_type", None))
        storage_temperature = attrs.get(
            "storage_temperature",
            getattr(self.instance, "storage_temperature", None),
        )

        item_type_name = (item_type.name if item_type else "").strip().lower()
        is_reagent = item_type_name in {"reagent", "reagents"}

        if is_reagent and not storage_temperature:
            raise serializers.ValidationError({
                "storage_temperature": "Storage temperature is required for reagents.",
            })

        return attrs

    def get_capacity_display(self, obj):
        if obj.capacity_value is None and not obj.capacity_unit:
            return None
        if obj.capacity_value is None:
            return obj.capacity_unit
        if obj.capacity_unit:
            return f"{obj.capacity_value} {obj.capacity_unit}"
        return str(obj.capacity_value)

    def get_label(self, obj):
        return obj.product_name

    def get_storage_temperature_label(self, obj):
        return obj.get_storage_temperature_display() or None
