from decimal import Decimal

from inventory.static_models import MaterialMaster


def infer_stock_unit_name(item_type_name, product_name):
    """
    Infer the stock unit from item type / product name.
    """
    item_type_name = (item_type_name or "").strip().lower()
    product_name_lower = (product_name or "").strip().lower()

    mapping = {
        "tip box": "box",
        "tube": "tube",
        "plate": "plate",
        "reservoir": "reservoir",
        "lid": "lid",
        "seal": "seal",
        "bottle": "bottle",
        "cap": "cap",
        "syringe": "syringe",
        "labels": "label",
        "reagent": "bottle",
        "serological pipette": "piece",
        "cell culture consumable": "piece",
        "spare part": "piece",
        "accessory": "piece",
        "other": "piece",
    }

    if item_type_name in mapping:
        return mapping[item_type_name]

    if "glove" in product_name_lower:
        return "box"

    if "mask" in product_name_lower:
        return "box"

    return "piece"


def infer_base_unit_name(item_type_name, product_name):
    """
    Infer the smallest logical unit for the material.
    """
    item_type_name = (item_type_name or "").strip().lower()
    product_name_lower = (product_name or "").strip().lower()

    mapping = {
        "tip box": "tip",
        "tube": "tube",
        "plate": "plate",
        "reservoir": "reservoir",
        "lid": "lid",
        "seal": "seal",
        "bottle": "bottle",
        "cap": "cap",
        "syringe": "syringe",
        "labels": "label",
        "reagent": "bottle",
        "serological pipette": "piece",
        "cell culture consumable": "piece",
        "spare part": "piece",
        "accessory": "piece",
        "other": "piece",
    }

    if item_type_name in mapping:
        return mapping[item_type_name]

    if "glove" in product_name_lower:
        return "glove"

    if "mask" in product_name_lower:
        return "mask"

    return "piece"


def calculate_base_units_per_stock_unit(single_pieces, units_per_box):
    """
    Excel packaging rule:
    Total Pieces = Single Pieces * Units per Box * Boxes in Stock

    Therefore:
    one stock unit contains Single Pieces * Units per Box base units
    """
    single_pieces = single_pieces or Decimal("1")
    units_per_box = units_per_box or Decimal("1")

    return single_pieces * units_per_box


def get_or_create_material_master(
    *,
    product_name,
    brand,
    manufacturer,
    vendor,
    manufacturer_catalog_number,
    vendor_catalog_number,
    item_type,
    device_type,
    capacity_value,
    capacity_unit,
):
    """
    Create or reuse a MaterialMaster.

    Identity is based on:
    - product_name
    - manufacturer_catalog_number
    - vendor_catalog_number
    """
    material, _ = MaterialMaster.objects.get_or_create(
        product_name=product_name,
        manufacturer_catalog_number=manufacturer_catalog_number,
        vendor_catalog_number=vendor_catalog_number,
        defaults={
            "brand": brand,
            "manufacturer": manufacturer,
            "vendor": vendor,
            "item_type": item_type,
            "device_type": device_type,
            "capacity_value": capacity_value,
            "capacity_unit": capacity_unit,
        },
    )

    fields_to_update = []

    if not material.brand and brand:
        material.brand = brand
        fields_to_update.append("brand")

    if not material.manufacturer and manufacturer:
        material.manufacturer = manufacturer
        fields_to_update.append("manufacturer")

    if not material.vendor and vendor:
        material.vendor = vendor
        fields_to_update.append("vendor")

    if not material.item_type and item_type:
        material.item_type = item_type
        fields_to_update.append("item_type")

    if not material.device_type and device_type:
        material.device_type = device_type
        fields_to_update.append("device_type")

    if material.capacity_value is None and capacity_value is not None:
        material.capacity_value = capacity_value
        fields_to_update.append("capacity_value")

    if not material.capacity_unit and capacity_unit:
        material.capacity_unit = capacity_unit
        fields_to_update.append("capacity_unit")

    if fields_to_update:
        material.save(update_fields=fields_to_update)

    return material