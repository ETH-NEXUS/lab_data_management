import csv
from decimal import Decimal
from pathlib import Path

from django.contrib.auth import get_user_model

from helpers.logger import logger
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
    to_decimal,
)


def _detect_csv_encoding(csv_path):
    """
    Detect the CSV encoding from a small byte sample.
    """
    raw_sample = csv_path.read_bytes()[:65536]

    for encoding in ("utf-8-sig", "utf-8", "cp1252", "latin-1"):
        try:
            raw_sample.decode(encoding)
            return encoding
        except UnicodeDecodeError:
            continue

    return "utf-8"


def _detect_csv_delimiter(sample_text):
    """
    Detect delimiter from text sample; fall back to a simple heuristic.
    """
    try:
        dialect = csv.Sniffer().sniff(sample_text, delimiters=",;\t|")
        return dialect.delimiter
    except csv.Error:
        if sample_text.count(";") > sample_text.count(","):
            return ";"
        return ","


def _normalize_row_keys(row):
    """
    Normalize row keys by stripping whitespace from column names.
    """
    normalized = {}

    for key, value in row.items():
        if isinstance(key, str):
            key = key.strip()
        normalized[key] = value

    return normalized


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


def _get_or_create_or_update_base_unit(material, base_uom, stock_unit_name, base_unit_name):
    """
    Ensure the material has a base unit row.
    """
    material_base_unit, created = MaterialUnit.objects.get_or_create(
        material=material,
        unit=base_uom,
        defaults={
            "is_base_unit": True,
            "is_stock_unit": stock_unit_name == base_unit_name,
            "is_order_unit": False,
            "base_units_per_unit": Decimal("1"),
        },
    )

    fields_to_update = []

    if not material_base_unit.is_base_unit:
        material_base_unit.is_base_unit = True
        fields_to_update.append("is_base_unit")

    expected_is_stock_unit = stock_unit_name == base_unit_name
    if material_base_unit.is_stock_unit != expected_is_stock_unit:
        material_base_unit.is_stock_unit = expected_is_stock_unit
        fields_to_update.append("is_stock_unit")

    if material_base_unit.base_units_per_unit != Decimal("1"):
        material_base_unit.base_units_per_unit = Decimal("1")
        fields_to_update.append("base_units_per_unit")

    if fields_to_update:
        material_base_unit.save(update_fields=fields_to_update)

    return material_base_unit


def _get_or_create_or_update_stock_unit(
    material,
    stock_uom,
    stock_unit_name,
    base_unit_name,
    base_units_per_stock_unit,
    single_pieces,
    units_per_box,
):
    """
    Ensure the material has a stock unit row and keep the conversion factor updated.
    """
    material_stock_unit, created = MaterialUnit.objects.get_or_create(
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

    fields_to_update = []

    expected_is_base_unit = stock_unit_name == base_unit_name
    if material_stock_unit.is_base_unit != expected_is_base_unit:
        material_stock_unit.is_base_unit = expected_is_base_unit
        fields_to_update.append("is_base_unit")

    if not material_stock_unit.is_stock_unit:
        material_stock_unit.is_stock_unit = True
        fields_to_update.append("is_stock_unit")

    if not material_stock_unit.is_order_unit:
        material_stock_unit.is_order_unit = True
        fields_to_update.append("is_order_unit")

    if material_stock_unit.base_units_per_unit != base_units_per_stock_unit:
        material_stock_unit.base_units_per_unit = base_units_per_stock_unit
        fields_to_update.append("base_units_per_unit")

    if not material_stock_unit.notes:
        material_stock_unit.notes = (
            f"Imported from Excel packaging. "
            f"single_pieces={single_pieces}, units_per_box={units_per_box}"
        )
        fields_to_update.append("notes")

    if fields_to_update:
        material_stock_unit.save(update_fields=fields_to_update)

    return material_stock_unit


def _upsert_inventory_stock(
    *,
    material,
    sector,
    stock_unit,
    quantity,
    minimum_quantity,
    lot_number,
    expiry_date,
    notes,
):
    """
    Create or update a stock row.

    This makes the import idempotent:
    running the same import again updates the same stock row
    instead of creating duplicates.
    """
    stock, created = InventoryStock.objects.update_or_create(
        material=material,
        sector=sector,
        stock_unit=stock_unit,
        lot_number=lot_number,
        expiry_date=expiry_date,
        defaults={
            "quantity": quantity,
            "minimum_quantity": minimum_quantity,
            "notes": notes,
        },
    )

    # Preserve old notes if the row already existed and new import is partial.
    if not created:
        updated_fields = []

        merged_notes = merge_notes(stock.notes, notes)
        if stock.notes != merged_notes:
            stock.notes = merged_notes
            updated_fields.append("notes")

        if updated_fields:
            stock.save(update_fields=updated_fields)

    return stock


def _find_user_by_username(username):
    """
    Resolve the order owner to a user if possible.
    """
    if not username:
        return None

    user_model = get_user_model()
    user = user_model.objects.filter(username=username).first()

    if not user:
        logger.warning(
            f"No user found for order owner '{username}'. Storing null user."
        )

    return user


def _upsert_order(
    *,
    material,
    order_unit,
    amount,
    order_date,
    status,
    who_ordered,
    project,
    notes,
):
    """
    Create or update an order row.

    Since the model currently has no dedicated external import ID,
    we use a stable business key that should remain identical
    across repeated imports of the same CSV.
    """
    order, created = Order.objects.update_or_create(
        material=material,
        order_unit=order_unit,
        amount=amount,
        order_date=order_date,
        status=status,
        who_ordered=who_ordered,
        project=project,
        defaults={
            "notes": notes,
        },
    )

    if not created:
        merged_notes = merge_notes(order.notes, notes)
        if order.notes != merged_notes:
            order.notes = merged_notes
            order.save(update_fields=["notes"])

    return order


def import_inventory_row(row):
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

    _get_or_create_or_update_base_unit(
        material=material,
        base_uom=base_uom,
        stock_unit_name=stock_unit_name,
        base_unit_name=base_unit_name,
    )

    base_units_per_stock_unit = calculate_base_units_per_stock_unit(
        single_pieces=single_pieces,
        units_per_box=units_per_box,
    )

    material_stock_unit = _get_or_create_or_update_stock_unit(
        material=material,
        stock_uom=stock_uom,
        stock_unit_name=stock_unit_name,
        base_unit_name=base_unit_name,
        base_units_per_stock_unit=base_units_per_stock_unit,
        single_pieces=single_pieces,
        units_per_box=units_per_box,
    )

    room = get_or_create_named(Room, room_name) if room_name else None
    sector = None

    if room and sector_name:
        sector, _ = Sector.objects.get_or_create(
            room=room,
            name=sector_name,
        )

    if sector:
        _upsert_inventory_stock(
            material=material,
            sector=sector,
            stock_unit=material_stock_unit,
            quantity=boxes_in_stock,
            minimum_quantity=minimum_amount,
            lot_number=lot_number,
            expiry_date=expiry_date,
            notes=notes,
        )

    if order_status_raw or who_ordered or assigned_project_name or order_date_raw:
        project = find_project_by_name(assigned_project_name)
        order_status = normalize_order_status(order_status_raw)
        who_ordered_user = _find_user_by_username(who_ordered)

        _upsert_order(
            material=material,
            order_unit=material_stock_unit,
            amount=boxes_in_stock if boxes_in_stock > 0 else Decimal("1"),
            order_date=order_date,
            status=order_status,
            who_ordered=who_ordered_user,
            project=project,
            notes=notes,
        )

    return True


def import_inventory_csv(csv_path):
    """
    Import a whole CSV file.

    Returns:
        dict with summary counts
    """
    csv_path = Path(csv_path)

    imported_rows = 0
    skipped_rows = 0

    encoding = _detect_csv_encoding(csv_path)

    with csv_path.open("r", encoding=encoding, newline="") as csv_file:
        sample = csv_file.read(8192)
        csv_file.seek(0)

        delimiter = _detect_csv_delimiter(sample)
        reader = csv.DictReader(csv_file, delimiter=delimiter)

        fieldnames = [name.strip() for name in (reader.fieldnames or []) if name]

        logger.info(
            f"Detected CSV format for {csv_path}: encoding={encoding}, delimiter='{delimiter}'."
        )

        if "Product Name" not in fieldnames:
            raise ValueError(
                "CSV header is missing required column 'Product Name'. "
                "Make sure you pass the inventory CSV file (not the JSON seed file)."
            )

        for row in reader:
            row = _normalize_row_keys(row)
            imported = import_inventory_row(row)

            if imported:
                imported_rows += 1
            else:
                skipped_rows += 1

    return {
        "imported_rows": imported_rows,
        "skipped_rows": skipped_rows,
    }
