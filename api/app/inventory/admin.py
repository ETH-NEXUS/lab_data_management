
from decimal import Decimal

from django.contrib import admin
from django.db.models import Count
from django.utils.html import format_html

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
from .dynamic_models import (
    InventoryStock,
    InventoryStockTablePreference,
    MaterialUsage,
    Order,
    Room,
    Sector,
)


# =========================================================
# Shared admin mixins / helpers
# =========================================================

class NameOnlyAdmin(admin.ModelAdmin):
    """
    Reusable admin for simple lookup tables with only a `name` field.
    """
    search_fields = ("name",)
    ordering = ("name",)
    list_display = ("name",)
    list_per_page = 50


# =========================================================
# Lookup / static reference models
# =========================================================

@admin.register(Manufacturer)
class ManufacturerAdmin(NameOnlyAdmin):
    list_display = ("name", "material_count")

    def get_queryset(self, request):
        queryset = super().get_queryset(request)
        return queryset.annotate(materials_count=Count("materials", distinct=True))

    @admin.display(description="Materials", ordering="materials_count")
    def material_count(self, obj):
        return obj.materials_count


@admin.register(Brand)
class BrandAdmin(NameOnlyAdmin):
    list_display = ("name", "material_count")

    def get_queryset(self, request):
        queryset = super().get_queryset(request)
        return queryset.annotate(materials_count=Count("materials", distinct=True))

    @admin.display(description="Materials", ordering="materials_count")
    def material_count(self, obj):
        return obj.materials_count


@admin.register(Vendor)
class VendorAdmin(NameOnlyAdmin):
    list_display = ("name", "material_count")

    def get_queryset(self, request):
        queryset = super().get_queryset(request)
        return queryset.annotate(materials_count=Count("materials", distinct=True))

    @admin.display(description="Materials", ordering="materials_count")
    def material_count(self, obj):
        return obj.materials_count


@admin.register(DeviceType)
class DeviceTypeAdmin(NameOnlyAdmin):
    list_display = ("name", "material_count")

    def get_queryset(self, request):
        queryset = super().get_queryset(request)
        return queryset.annotate(materials_count=Count("materials", distinct=True))

    @admin.display(description="Materials", ordering="materials_count")
    def material_count(self, obj):
        return obj.materials_count


@admin.register(ItemType)
class ItemTypeAdmin(NameOnlyAdmin):
    list_display = ("name", "material_count")

    def get_queryset(self, request):
        queryset = super().get_queryset(request)
        return queryset.annotate(materials_count=Count("materials", distinct=True))

    @admin.display(description="Materials", ordering="materials_count")
    def material_count(self, obj):
        return obj.materials_count


@admin.register(MaterialAttribute)
class MaterialAttributeAdmin(NameOnlyAdmin):
    list_display = ("name", "material_count")

    def get_queryset(self, request):
        queryset = super().get_queryset(request)
        return queryset.annotate(materials_count=Count("materials", distinct=True))

    @admin.display(description="Materials", ordering="materials_count")
    def material_count(self, obj):
        return obj.materials_count


@admin.register(UnitOfMeasure)
class UnitOfMeasureAdmin(NameOnlyAdmin):
    list_display = ("name", "material_unit_count")

    def get_queryset(self, request):
        queryset = super().get_queryset(request)
        return queryset.annotate(material_units_count=Count("material_units", distinct=True))

    @admin.display(description="Material units", ordering="material_units_count")
    def material_unit_count(self, obj):
        return obj.material_units_count


# =========================================================
# Inlines
# =========================================================

class MaterialUnitInline(admin.TabularInline):
    model = MaterialUnit
    extra = 1
    autocomplete_fields = ("unit",)
    fields = (
        "unit",
        "is_base_unit",
        "is_stock_unit",
        "is_order_unit",
        "base_units_per_unit",
        "notes",
    )
    verbose_name = "Material unit"
    verbose_name_plural = "Material units"


class InventoryStockInline(admin.TabularInline):
    model = InventoryStock
    extra = 0
    autocomplete_fields = ("sector", "stock_unit")
    fields = (
        "sector",
        "stock_unit",
        "quantity",
        "minimum_quantity",
        "lot_number",
        "expiry_date",
        "is_favorite",
        "notes",
    )
    verbose_name = "Stock entry"
    verbose_name_plural = "Stock entries"
    show_change_link = True


class OrderInline(admin.TabularInline):
    model = Order
    extra = 0
    autocomplete_fields = ("order_unit", "project", "who_ordered")
    fields = (
        "order_date",
        "status",
        "order_unit",
        "amount",
        "who_ordered",
        "project",
        "notes",
    )
    verbose_name = "Order"
    verbose_name_plural = "Orders"
    show_change_link = True


class MaterialUsageInline(admin.TabularInline):
    model = MaterialUsage
    extra = 0
    autocomplete_fields = ("project", "experiment", "usage_unit")
    fields = (
        "project",
        "experiment",
        "usage_unit",
        "quantity_used",
        "used_at",
        "notes",
    )
    readonly_fields = ("used_at",)
    verbose_name = "Material usage"
    verbose_name_plural = "Material usages"
    show_change_link = True


class SectorInline(admin.TabularInline):
    model = Sector
    extra = 1
    fields = ("name",)
    verbose_name = "Sector"
    verbose_name_plural = "Sectors"


# =========================================================
# Material master admin
# =========================================================

@admin.action(description="Mark selected materials as active")
def mark_materials_active(modeladmin, request, queryset):
    queryset.update(is_active=True)


@admin.action(description="Mark selected materials as inactive")
def mark_materials_inactive(modeladmin, request, queryset):
    queryset.update(is_active=False)


@admin.register(MaterialMaster)
class MaterialMasterAdmin(admin.ModelAdmin):
    list_display = (
        "product_name",
        "brand",
        "manufacturer",
        "vendor",
        "item_type",
        "device_type",
        "capacity_display",
        "attribute_list",
        "is_active",
        "stock_entry_count",
    )
    list_filter = (
        "is_active",
        "item_type",
        "device_type",
        "brand",
        "manufacturer",
        "vendor",
        "attributes",
    )
    search_fields = (
        "product_name",
        "manufacturer_catalog_number",
        "vendor_catalog_number",
        "serial_number",
        "order_number",
        "description",
        "brand__name",
        "manufacturer__name",
        "vendor__name",
        "item_type__name",
        "device_type__name",
        "attributes__name",
    )
    autocomplete_fields = (
        "brand",
        "manufacturer",
        "vendor",
        "item_type",
        "device_type",
    )
    filter_horizontal = ("attributes",)
    readonly_fields = (
        "stock_entry_count",
        "total_stock_base_units",
        "unit_summary",
    )
    list_per_page = 50
    actions = (mark_materials_active, mark_materials_inactive)
    inlines = (MaterialUnitInline, InventoryStockInline, OrderInline)

    fieldsets = (
        (
            "Core product information",
            {
                "fields": (
                    "product_name",
                    "brand",
                    "manufacturer",
                    "vendor",
                    "item_type",
                    "device_type",
                    "attributes",
                    "is_active",
                )
            },
        ),
        (
            "Catalog references",
            {
                "fields": (
                    "manufacturer_catalog_number",
                    "vendor_catalog_number",
                    "serial_number",
                    "order_number",
                )
            },
        ),
        (
            "Technical details",
            {
                "fields": (
                    "capacity_value",
                    "capacity_unit",
                    "storage_temperature",
                    "safety_data_sheet",
                    "default_cost",
                    "lifetime_days",
                    "description",
                )
            },
        ),
        (
            "Quick inventory summary",
            {
                "classes": ("collapse",),
                "fields": (
                    "stock_entry_count",
                    "total_stock_base_units",
                    "unit_summary",
                ),
            },
        ),
    )

    def get_queryset(self, request):
        queryset = super().get_queryset(request)
        return queryset.annotate(stock_entries_count=Count("stock_entries", distinct=True))

    @admin.display(description="Capacity")
    def capacity_display(self, obj):
        if obj.capacity_value is None and not obj.capacity_unit:
            return "-"
        if obj.capacity_value is None:
            return obj.capacity_unit
        if obj.capacity_unit:
            return f"{obj.capacity_value} {obj.capacity_unit}"
        return str(obj.capacity_value)

    @admin.display(description="Attributes")
    def attribute_list(self, obj):
        attributes = obj.attributes.all()
        if not attributes:
            return "-"
        return ", ".join(attribute.name for attribute in attributes)

    @admin.display(description="Stock entries", ordering="stock_entries_count")
    def stock_entry_count(self, obj):
        return obj.stock_entries_count

    @admin.display(description="Total stock in base units")
    def total_stock_base_units(self, obj):
        total = Decimal("0")
        for stock in obj.stock_entries.select_related("stock_unit").all():
            total += stock.quantity_in_base_units
        return total

    @admin.display(description="Units summary")
    def unit_summary(self, obj):
        units = obj.units.select_related("unit").all()
        if not units:
            return "-"
        return "; ".join(
            f"{unit.unit.name}: {unit.base_units_per_unit}"
            for unit in units
        )


# =========================================================
# Material unit admin
# =========================================================

@admin.register(MaterialUnit)
class MaterialUnitAdmin(admin.ModelAdmin):
    list_display = (
        "material",
        "unit",
        "is_base_unit",
        "is_stock_unit",
        "is_order_unit",
        "base_units_per_unit",
    )
    list_filter = (
        "is_base_unit",
        "is_stock_unit",
        "is_order_unit",
        "unit",
        "material__item_type",
    )
    search_fields = (
        "material__product_name",
        "unit__name",
        "notes",
    )
    autocomplete_fields = ("material", "unit")
    list_editable = (
        "is_base_unit",
        "is_stock_unit",
        "is_order_unit",
        "base_units_per_unit",
    )
    list_per_page = 100


# =========================================================
# Location admin
# =========================================================

@admin.register(Room)
class RoomAdmin(admin.ModelAdmin):
    list_display = ("name", "sector_count", "stock_entry_count")
    search_fields = ("name",)
    inlines = (SectorInline,)
    list_per_page = 50

    def get_queryset(self, request):
        queryset = super().get_queryset(request)
        return queryset.annotate(
            sectors_count=Count("sectors", distinct=True),
            stock_entries_count=Count("sectors__stock_entries", distinct=True),
        )

    @admin.display(description="Sectors", ordering="sectors_count")
    def sector_count(self, obj):
        return obj.sectors_count

    @admin.display(description="Stock entries", ordering="stock_entries_count")
    def stock_entry_count(self, obj):
        return obj.stock_entries_count


@admin.register(Sector)
class SectorAdmin(admin.ModelAdmin):
    list_display = ("name", "room", "stock_entry_count")
    list_filter = ("room",)
    search_fields = ("name", "room__name")
    autocomplete_fields = ("room",)
    list_per_page = 100

    def get_queryset(self, request):
        queryset = super().get_queryset(request)
        return queryset.annotate(stock_entries_count=Count("stock_entries", distinct=True))

    @admin.display(description="Stock entries", ordering="stock_entries_count")
    def stock_entry_count(self, obj):
        return obj.stock_entries_count


# =========================================================
# Inventory stock admin
# =========================================================

@admin.action(description="Mark selected stock entries as favorite")
def mark_stock_favorite(modeladmin, request, queryset):
    queryset.update(is_favorite=True)


@admin.action(description="Unmark selected stock entries as favorite")
def unmark_stock_favorite(modeladmin, request, queryset):
    queryset.update(is_favorite=False)


@admin.register(InventoryStock)
class InventoryStockAdmin(admin.ModelAdmin):
    list_display = (
        "material",
        "room_display",
        "sector",
        "stock_unit",
        "quantity",
        "minimum_quantity",
        "inventory_status_colored",
        "quantity_in_base_units_display",
        "lot_number",
        "expiry_date",
        "is_favorite",
        "updated_at",
    )
    list_filter = (
        "is_favorite",
        "sector__room",
        "material__item_type",
        "material__device_type",
        "material__manufacturer",
        "material__vendor",
        "expiry_date",
        "created_at",
        "updated_at",
    )
    search_fields = (
        "material__product_name",
        "material__manufacturer_catalog_number",
        "material__vendor_catalog_number",
        "lot_number",
        "notes",
        "sector__name",
        "sector__room__name",
    )
    autocomplete_fields = ("material", "sector", "stock_unit")
    readonly_fields = (
        "created_at",
        "updated_at",
        "inventory_status",
        "quantity_in_base_units_display",
        "minimum_quantity_in_base_units_display",
    )
    list_editable = (
        "quantity",
        "minimum_quantity",
        "is_favorite",
    )
    actions = (mark_stock_favorite, unmark_stock_favorite)
    list_per_page = 100
    date_hierarchy = "updated_at"

    fieldsets = (
        (
            "Material and location",
            {
                "fields": (
                    "material",
                    "sector",
                    "stock_unit",
                )
            },
        ),
        (
            "Stock amounts",
            {
                "fields": (
                    "quantity",
                    "minimum_quantity",
                    "inventory_status",
                    "quantity_in_base_units_display",
                    "minimum_quantity_in_base_units_display",
                )
            },
        ),
        (
            "Batch information",
            {
                "fields": (
                    "lot_number",
                    "expiry_date",
                )
            },
        ),
        (
            "Extra",
            {
                "fields": (
                    "is_favorite",
                    "notes",
                    "created_at",
                    "updated_at",
                )
            },
        ),
    )

    inlines = (MaterialUsageInline,)

    @admin.display(description="Room", ordering="sector__room__name")
    def room_display(self, obj):
        return obj.sector.room if obj.sector_id else "-"

    @admin.display(description="Status")
    def inventory_status_colored(self, obj):
        if obj.inventory_status == "low":
            return format_html('<b style="color:#b91c1c;">Low</b>')
        return format_html('<span style="color:#15803d;">In stock</span>')

    @admin.display(description="Qty in base units")
    def quantity_in_base_units_display(self, obj):
        return obj.quantity_in_base_units

    @admin.display(description="Min in base units")
    def minimum_quantity_in_base_units_display(self, obj):
        return obj.minimum_quantity_in_base_units


# =========================================================
# Order admin
# =========================================================

@admin.register(Order)
class OrderAdmin(admin.ModelAdmin):
    list_display = (
        "material",
        "amount",
        "order_unit",
        "status",
        "who_ordered",
        "project",
        "order_date",
        "created_at",
    )
    list_filter = (
        "status",
        "order_unit",
        "material__item_type",
        "material__manufacturer",
        "material__vendor",
        "order_date",
        "created_at",
    )
    search_fields = (
        "material__product_name",
        "material__manufacturer_catalog_number",
        "material__vendor_catalog_number",
        "notes",
        "project__name",
        "who_ordered__username",
        "who_ordered__first_name",
        "who_ordered__last_name",
        "who_ordered__email",
    )
    autocomplete_fields = (
        "material",
        "order_unit",
        "project",
        "who_ordered",
    )
    readonly_fields = ("created_at", "updated_at")
    list_editable = ("status",)
    list_per_page = 100
    date_hierarchy = "order_date"

    fieldsets = (
        (
            "Order details",
            {
                "fields": (
                    "material",
                    "order_unit",
                    "amount",
                    "status",
                    "order_date",
                )
            },
        ),
        (
            "Who / project",
            {
                "fields": (
                    "who_ordered",
                    "project",
                )
            },
        ),
        (
            "Notes and timestamps",
            {
                "fields": (
                    "notes",
                    "created_at",
                    "updated_at",
                )
            },
        ),
    )


# =========================================================
# Material usage admin
# =========================================================

@admin.register(MaterialUsage)
class MaterialUsageAdmin(admin.ModelAdmin):
    list_display = (
        "project",
        "experiment",
        "material_display",
        "inventory_stock",
        "usage_unit",
        "quantity_used",
        "quantity_used_in_base_units_display",
        "used_at",
    )
    list_filter = (
        "project",
        "experiment",
        "usage_unit",
        "inventory_stock__material__item_type",
        "inventory_stock__material__device_type",
        "used_at",
    )
    search_fields = (
        "project__name",
        "experiment__name",
        "inventory_stock__material__product_name",
        "inventory_stock__lot_number",
        "notes",
    )
    autocomplete_fields = (
        "project",
        "experiment",
        "inventory_stock",
        "usage_unit",
    )
    readonly_fields = (
        "used_at",
        "quantity_used_in_base_units_display",
    )
    list_per_page = 100
    date_hierarchy = "used_at"

    fieldsets = (
        (
            "Usage record",
            {
                "fields": (
                    "project",
                    "experiment",
                    "inventory_stock",
                    "usage_unit",
                    "quantity_used",
                    "quantity_used_in_base_units_display",
                    "used_at",
                )
            },
        ),
        (
            "Notes",
            {
                "fields": ("notes",)
            },
        ),
    )

    @admin.display(description="Material")
    def material_display(self, obj):
        if obj.inventory_stock_id:
            return obj.inventory_stock.material
        return "-"

    @admin.display(description="Qty in base units")
    def quantity_used_in_base_units_display(self, obj):
        return obj.quantity_used_in_base_units


# =========================================================
# Table preference admin
# =========================================================

@admin.register(InventoryStockTablePreference)
class InventoryStockTablePreferenceAdmin(admin.ModelAdmin):
    list_display = (
        "user",
        "table_key",
        "preset",
        "updated_at",
    )
    list_filter = (
        "table_key",
        "preset",
        "updated_at",
    )
    search_fields = (
        "user__username",
        "user__first_name",
        "user__last_name",
        "user__email",
    )
    autocomplete_fields = ("user",)
    readonly_fields = (
        "created_at",
        "updated_at",
    )
    list_per_page = 100
