from django.db import models
from django.db.models import Q

from .dynamic_models import InventoryStock


def build_inventory_stock_base_queryset():
    """
    Returns the shared stock queryset used by paginated inventory list endpoints.

    Returned queryset example:
    - `InventoryStock.objects.select_related(...).prefetch_related(...)`
    """
    return (
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
            "source_order",
        )
        .prefetch_related(
            "material__attributes",
            "material__units__unit",
            "additional_sectors",
            "additional_sectors__room",
        )
    )


def apply_inventory_stock_list_filters(queryset, query_params, favorite_user=None):
    """
    Applies the shared stock-table filters and ordering for stock list endpoints.

    Accepted query param examples:
    - `{ "search": "pbs", "ordering": "-created_at" }`
    - `{ "room": "2", "inventory_status": "low" }`

    Returned queryset example:
    - `queryset.filter(...).annotate(...).order_by("-created_at")`
    """
    filtered_queryset = queryset.annotate(
        inventory_status_rank=models.Case(
            models.When(quantity=0, then=models.Value(2)),
            models.When(quantity__lt=models.F("minimum_quantity"), then=models.Value(1)),
            default=models.Value(0),
            output_field=models.IntegerField(),
        ),
    )

    search = query_params.get("search")
    room_id = query_params.get("room")
    sector_id = query_params.get("sector")
    item_type_id = query_params.get("item_type")
    device_type_id = query_params.get("device_type")
    manufacturer_id = query_params.get("manufacturer")
    vendor_id = query_params.get("vendor")
    is_favorite = query_params.get("is_favorite")
    inventory_status = query_params.get("inventory_status")
    lot_number = query_params.get("lot_number")
    ordering = query_params.get("ordering")

    if search:
        filtered_queryset = filtered_queryset.filter(
            Q(material__product_name__icontains=search)
            | Q(material__manufacturer_catalog_number__icontains=search)
            | Q(material__vendor_catalog_number__icontains=search)
            | Q(lot_number__icontains=search)
            | Q(notes__icontains=search)
            | Q(sector__name__icontains=search)
            | Q(sector__room__name__icontains=search)
            | Q(additional_sectors__name__icontains=search)
            | Q(additional_sectors__room__name__icontains=search)
        ).distinct()

    if room_id:
        filtered_queryset = filtered_queryset.filter(
            Q(sector__room_id=room_id) | Q(additional_sectors__room_id=room_id)
        ).distinct()

    if sector_id:
        filtered_queryset = filtered_queryset.filter(
            Q(sector_id=sector_id) | Q(additional_sectors__id=sector_id)
        ).distinct()

    if item_type_id:
        filtered_queryset = filtered_queryset.filter(material__item_type_id=item_type_id)

    if device_type_id:
        filtered_queryset = filtered_queryset.filter(material__device_type_id=device_type_id)

    if manufacturer_id:
        filtered_queryset = filtered_queryset.filter(material__manufacturer_id=manufacturer_id)

    if vendor_id:
        filtered_queryset = filtered_queryset.filter(material__vendor_id=vendor_id)

    if is_favorite is not None:
        wants_favorites = is_favorite.lower() in ("true", "1", "yes")

        if favorite_user is None or not favorite_user.is_authenticated:
            if wants_favorites:
                return filtered_queryset.none()
        elif wants_favorites:
            filtered_queryset = filtered_queryset.filter(favorite_users=favorite_user)
        else:
            filtered_queryset = filtered_queryset.exclude(favorite_users=favorite_user)

    if inventory_status == "low":
        filtered_queryset = filtered_queryset.filter(
            quantity__gt=0,
            quantity__lt=models.F("minimum_quantity"),
        )

    if inventory_status == "out_of_stock":
        filtered_queryset = filtered_queryset.filter(quantity=0)

    if inventory_status == "in_stock":
        filtered_queryset = filtered_queryset.filter(quantity__gte=models.F("minimum_quantity"))

    if lot_number:
        filtered_queryset = filtered_queryset.filter(lot_number__icontains=lot_number)

    if ordering:
        ordering_fields = []
        ordering_map = {
            "productName": ["material__product_name"],
            "favorite": ["is_favorite"],
            "inventoryStatus": ["inventory_status_rank"],
            "quantityWithStockUnit": ["quantity"],
            "minimumQuantity": ["minimum_quantity"],
            "location": ["sector__room__name", "sector__name"],
            "deviceType": ["material__device_type__name"],
            "itemType": ["material__item_type__name"],
            "storageTemperature": ["material__storage_temperature"],
            "brand": ["material__brand__name"],
            "manufacturer": ["material__manufacturer__name"],
            "vendor": ["material__vendor__name"],
            "manufacturerCatalogNumber": ["material__manufacturer_catalog_number"],
            "vendorCatalogNumber": ["material__vendor_catalog_number"],
            "capacity": ["material__capacity_value"],
            "defaultCost": ["material__default_cost"],
            "isActive": ["material__is_active"],
            "description": ["material__description"],
            "serialNumber": ["material__serial_number"],
            "orderNumber": ["material__order_number"],
            "lifetimeDays": ["material__lifetime_days"],
            "lotNumber": ["lot_number"],
            "expiryDate": ["expiry_date"],
            "archivedAt": ["archived_at"],
            "notes": ["notes"],
            "createdAt": ["created_at"],
            "updatedAt": ["updated_at"],
        }

        for raw_value in ordering.split(","):
            normalized_value = raw_value.strip()

            if normalized_value == "":
                continue

            is_descending = normalized_value.startswith("-")
            ordering_key = normalized_value[1:] if is_descending else normalized_value
            mapped_fields = ordering_map.get(ordering_key)

            if not mapped_fields:
                continue

            for mapped_field in mapped_fields:
                ordering_fields.append(f"-{mapped_field}" if is_descending else mapped_field)

        if ordering_fields:
            return filtered_queryset.order_by(*ordering_fields)

    return filtered_queryset.order_by("-created_at")
