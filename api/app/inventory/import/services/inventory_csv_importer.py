import csv
from decimal import Decimal
from pathlib import Path

from inventory.dynamic_models import InventoryStock, Order, Room, Sector
from inventory.static_models import (
    Brand,
    DeviceType,
    ItemType,
    Manufacturer,
    MaterialAttribute,
    MaterialUnit,
    UnitOfMeasure,
    Vendor,
)
from ..utils.dates import parse_flexible_date, parse_flexible_datetime
from ..utils.lookups import (
    find_project_by_name,
    get_or_create_named,
    normalize_order_status,
)
from ..utils.materials import (
    calculate_base_units_per_stock_unit,
    get_or_create_material_master,
    infer_base_unit_name,
    infer_stock_unit_name,
)
from ..utils.parsing import (
    clean_value,
    merge_notes,
    split_attributes,
    to_bool,
    to_decimal,
)


def should_skip_inventory_row(row):
    """
    Decide whether a CSV row should be skipped.
    """
    product_name = clean_value(row.get("Product Name"))

    if not product_name:
        return True

    lower_name = product_name.lower()

    helper_phrases = [
        "copy name from vendor",
        "product information",
        "stock information",
        "order information",
    ]

    return any(phrase in lower_name for phrase in helper_phrases)


def import_inventory_row(row, *, update_existing_stock=False):
    """
    Import one CSV row into the database.

    Returns:
        True if imported
        False if skipped
    """
    if should_skip_inventory_row(row):
        return False

    product_name = clean_value(row.get("Product Name"))
    brand_name = clean_value(row.get("Brand"))
    manufacturer_name = clean_value(row.get("Manufacturer"))
    vendor_name = clean_value(row.get("Vendor"))
    manufacturer_cat_no = clean_value(row.get("Cat# Manufacturer"))
    vendor_cat_no = clean_value(row.get("Cat# Vendor"))
    lot_number = clean_value(row.get("Lot#"))
    expiry_date_raw = clean_value(row.get("Expiry Date"))
    room_name = clean_value(row.get("Room"))
    sector_name = clean_value(row.get("Sector"))
    capacity_unit = clean_value(row.get("Unit"))
    capacity_value = to_decimal(row.get("Capacity"))
    item_property_text = clean_value(row.get("Item Property"))
    device_name = clean_value(row.get("Pipetting Device"))
    item_type_name = clean_value(row.get("Item Type"))
    single_pieces = to_decimal(row.get("Single Pieces"), default=Decimal("1"))
    units_per_box = to_decimal(row.get("Units per Box"), default=Decimal("1"))
    boxes_in_stock = to_decimal(row.get("Boxes in Stock"), default=Decimal("0"))
    minimum_amount = to_decimal(row.get("Minimum Amount"), default=Decimal("0"))
    order_status_raw = clean_value(row.get("Order Status"))
    who_ordered = clean_value(row.get("Who Ordered"))
    assigned_project_name = clean_value(row.get("Assigned Project"))
    order_date_raw = clean_value(row.get("Order Date"))
    notes = clean_value(row.get("Notes"))
    favorites_raw = clean_value(row.get("Favorites"))

    expiry_date = parse_flexible_date(expiry_date_raw)
    order_date = parse_flexible_datetime(order_date_raw)

    brand = get_or_create_named(Brand, brand_name)
    manufacturer = get_or_create_named(Manufacturer, manufacturer_name)
    vendor = get_or_create_named(Vendor, vendor_name)
    device_type = get_or_create_named(DeviceType, device_name)
    item_type = get_or_create_named(ItemType, item_type_name)

    material = get_or_create_material_master(
        product_name=product_name,
        brand=brand,
        manufacturer=manufacturer,
        vendor=vendor,
        manufacturer_catalog_number=manufacturer_cat_no,
        vendor_catalog_number=vendor_cat_no,
        item_type=item_type,
        device_type=device_type,
        capacity_value=capacity_value,
        capacity_unit=capacity_unit,
    )

    for attribute_name in split_attributes(item_property_text):
        attribute = get_or_create_named(MaterialAttribute, attribute_name)
        if attribute:
            material.attributes.add(attribute)

    stock_unit_name = infer_stock_unit_name(item_type_name, product_name)
    base_unit_name = infer_base_unit_name(item_type_name, product_name)

    stock_uom = get_or_create_named(UnitOfMeasure, stock_unit_name)
    base_uom = get_or_create_named(UnitOfMeasure, base_unit_name)

    MaterialUnit.objects.get_or_create(
        material=material,
        unit=base_uom,
        defaults={
            "is_base_unit": True,
            "is_stock_unit": stock_unit_name == base_unit_name,
            "is_order_unit": False,
            "base_units_per_unit": Decimal("1"),
        },
    )

    base_units_per_stock_unit = calculate_base_units_per_stock_unit(
        single_pieces=single_pieces,
        units_per_box=units_per_box,
    )

    material_stock_unit, _ = MaterialUnit.objects.get_or_create(
        material=material,
        unit=stock_uom,
        defaults={
            "is_base_unit": stock_unit_name == base_unit_name,
            "is_stock_unit": True,
            "is_order_unit": True,
            "base_units_per_unit": base_units_per_stock_unit,
            "notes": (
                f"Imported from Excel packaging. "
                f"single_pieces={single_pieces}, units_per_box={units_per_box}"
            ),
        },
    )

    if material_stock_unit.base_units_per_unit in [None, Decimal("0")]:
        material_stock_unit.base_units_per_unit = base_units_per_stock_unit
        material_stock_unit.save(update_fields=["base_units_per_unit"])

    room = get_or_create_named(Room, room_name) if room_name else None
    sector = None

    if room and sector_name:
        sector, _ = Sector.objects.get_or_create(
            room=room,
            name=sector_name,
        )

    if sector:
        is_favorite = to_bool(favorites_raw)

        if update_existing_stock:
            stock = InventoryStock.objects.filter(
                material=material,
                sector=sector,
                stock_unit=material_stock_unit,
                lot_number=lot_number,
                expiry_date=expiry_date,
            ).first()

            if stock:
                stock.quantity = boxes_in_stock
                stock.minimum_quantity = minimum_amount
                stock.notes = merge_notes(stock.notes, notes)
                stock.is_favorite = stock.is_favorite or is_favorite
                stock.save()
            else:
                InventoryStock.objects.create(
                    material=material,
                    sector=sector,
                    stock_unit=material_stock_unit,
                    quantity=boxes_in_stock,
                    minimum_quantity=minimum_amount,
                    lot_number=lot_number,
                    expiry_date=expiry_date,
                    notes=notes,
                    is_favorite=is_favorite,
                )
        else:
            InventoryStock.objects.create(
                material=material,
                sector=sector,
                stock_unit=material_stock_unit,
                quantity=boxes_in_stock,
                minimum_quantity=minimum_amount,
                lot_number=lot_number,
                expiry_date=expiry_date,
                notes=notes,
                is_favorite=is_favorite,
            )

    if order_status_raw or who_ordered or assigned_project_name or order_date_raw:
        project = find_project_by_name(assigned_project_name)
        order_status = normalize_order_status(order_status_raw)

        Order.objects.create(
            material=material,
            order_unit=material_stock_unit,
            amount=boxes_in_stock if boxes_in_stock > 0 else Decimal("1"),
            order_date=order_date,
            status=order_status,
            who_ordered=who_ordered,
            project=project,
            notes=notes,
        )

    return True


def import_inventory_csv(csv_path, *, update_existing_stock=False):
    """
    Import a whole CSV file.

    Returns:
        dict with summary counts
    """
    csv_path = Path(csv_path)

    imported_rows = 0
    skipped_rows = 0

    with csv_path.open("r", encoding="utf-8-sig", newline="") as csv_file:
        reader = csv.DictReader(csv_file)

        for row in reader:
            imported = import_inventory_row(
                row,
                update_existing_stock=update_existing_stock,
            )

            if imported:
                imported_rows += 1
            else:
                skipped_rows += 1

    return {
        "imported_rows": imported_rows,
        "skipped_rows": skipped_rows,
    }
