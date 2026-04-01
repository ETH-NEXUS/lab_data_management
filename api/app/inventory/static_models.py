from django.db import models


class Manufacturer(models.Model):
    """
    Company that manufactures the product.

    - Tecan
    - Integra
    - Eppendorf
    - Corning
    - Greiner bio-one
    - Beckman Coulter
    - HighRes Biosolutions
    """
    name = models.CharField(
        max_length=255,
        unique=True,
        help_text="Manufacturer name, e.g. Tecan, Eppendorf, Corning.",
    )

    class Meta:
        ordering = ["name"]
        verbose_name = "Manufacturer"
        verbose_name_plural = "Manufacturers"

    def __str__(self):
        return self.name


class Brand(models.Model):
    """
    Optional brand name.

    - Costar
    - Cellstar
    - Twin.tec
    - fisherbrand
    - Sempercare
    - oxo
    """
    name = models.CharField(
        max_length=255,
        unique=True,
        help_text="Brand name if applicable, e.g. Costar, Cellstar, Twin.tec.",
    )

    class Meta:
        ordering = ["name"]
        verbose_name = "Brand"
        verbose_name_plural = "Brands"

    def __str__(self):
        return self.name


class Vendor(models.Model):
    """
    Supplier / seller from whom the lab orders the material.

    - Tecan
    - Witec
    - D-BIOL shop
    - Huberlab
    - Brack
    - Galaxus
    """
    name = models.CharField(
        max_length=255,
        unique=True,
        help_text="Vendor / supplier name, e.g. Witec, D-BIOL shop, Huberlab.",
    )

    class Meta:
        ordering = ["name"]
        verbose_name = "Vendor"
        verbose_name_plural = "Vendors"

    def __str__(self):
        return self.name


class DeviceType(models.Model):
    """
    Compatible pipetting or automation device.

    Examples from the lab's excel:
    - LiHa
    - MCA96
    - MCA384
    - Hand pipette

    Must be nullable on MaterialMaster because many items
    do not belong to a specific device.
    """
    name = models.CharField(
        max_length=255,
        unique=True,
        help_text="Compatible device type, e.g. LiHa, MCA96, MCA384, Hand pipette.",
    )

    class Meta:
        ordering = ["name"]
        verbose_name = "Device type"
        verbose_name_plural = "Device types"

    def __str__(self):
        return self.name


class ItemType(models.Model):
    """
    General category of the material.

    Examples from the lab's excel:
    - tip box
    - tube
    - plate
    - reservoir
    - lid
    - seal
    - bottle
    - syringe
    - reagent
    - spare part
    - accessory
    - cell culture consumable
    - serological pipette
    - cap
    - labels
    - other
    """
    name = models.CharField(
        max_length=255,
        unique=True,
        help_text="General material category, e.g. tip box, plate, tube, reagent, spare part.",
    )

    class Meta:
        ordering = ["name"]
        verbose_name = "Item type"
        verbose_name_plural = "Item types"

    def __str__(self):
        return self.name


class MaterialAttribute(models.Model):
    """
    Reusable product attribute.

    This is intentionally generic because the lab's excel file contains many
    repeatable descriptors that do not fit into one rigid schema.

    Examples from the lab's excel:
    - sterile
    - filtered
    - pure
    - conductive
    - wide bore
    - low retention
    - PCR clean
    - refill
    - single wrapped
    - individually wrapped
    - culture treated
    - non-treated
    - barcoded
    - RNase/DNase/protease/pyrogen free
    - non-filtered
    - non-sterile
    """
    name = models.CharField(
        max_length=255,
        unique=True,
        help_text="Reusable attribute, e.g. sterile, filtered, low retention, wide bore.",
    )

    class Meta:
        ordering = ["name"]
        verbose_name = "Material attribute"
        verbose_name_plural = "Material attributes"

    def __str__(self):
        return self.name


class UnitOfMeasure(models.Model):
    """
    Unit in which a material can be counted.

    This is for inventory counting / packaging logic,
    not for physical capacity.

    Examples from the lab's excel:
    - piece
    - box
    - bag
    - rack
    - blister
    - bottle
    - tip
    - tube
    - plate
    - lid
    - roll
    - cassette
    - pack
    """
    name = models.CharField(
        max_length=100,
        unique=True,
        help_text="Inventory unit, e.g. piece, box, bag, rack, plate, bottle.",
    )

    class Meta:
        ordering = ["name"]
        verbose_name = "Unit of measure"
        verbose_name_plural = "Units of measure"

    def __str__(self):
        return self.name


class MaterialMaster(models.Model):
    """
    One row = one real catalog product / stock keeping unit / material definition.

    This is master data and should describe the product itself,
    not a specific stock batch or shelf record.

    Examples from the lab's excel:
    - 200ul MCA96 Tips, pure, no filter
    - 1000ul LiHa, conductive, pure, wide bore, filter, 24 racks in 2er Blister
    - 50mL Centrifuge Tubes
    - Costar Assay Plate 96 Well, flat bottom, culture treated, with lid, barcoded, sterile
    - Gloves M, Nitril, Sempercare Safe
    - Prime spare power cable

    Things that should not live here:
    - lot number
    - expiry date of one batch
    - room / sector
    - boxes currently in stock
    - who ordered one shipment
    """

    product_name = models.CharField(
        max_length=255,
        help_text="Catalog product name copied from vendor/manufacturer if possible.",
    )

    brand = models.ForeignKey(
        Brand,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="materials",
        help_text="Optional brand, e.g. Costar, Cellstar, Twin.tec.",
    )

    manufacturer = models.ForeignKey(
        Manufacturer,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="materials",
        help_text="Manufacturer of the product.",
    )

    vendor = models.ForeignKey(
        Vendor,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="materials",
        help_text="Vendor / supplier from whom the item is ordered.",
    )

    manufacturer_catalog_number = models.CharField(
        max_length=255,
        null=True,
        blank=True,
        db_index=True,
        help_text="Manufacturer catalog number, e.g. 30038616, 3815, 710180.",
    )

    vendor_catalog_number = models.CharField(
        max_length=255,
        null=True,
        blank=True,
        db_index=True,
        help_text="Vendor catalog number. May be identical to manufacturer catalog number.",
    )

    item_type = models.ForeignKey(
        ItemType,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="materials",
        help_text="General item category, e.g. tip box, plate, tube, reagent, spare part.",
    )

    device_type = models.ForeignKey(
        DeviceType,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="materials",
        help_text="Compatible device if applicable, e.g. LiHa, MCA96, MCA384, Hand pipette.",
    )

    attributes = models.ManyToManyField(
        MaterialAttribute,
        blank=True,
        related_name="materials",
        help_text="Flexible repeatable attributes like sterile, filtered, pure, wide bore.",
    )

    capacity_value = models.DecimalField(
        max_digits=12,
        decimal_places=3,
        null=True,
        blank=True,
        help_text="Capacity / size of one single piece, e.g. 200, 50, 1.5, 75.",
    )

    capacity_unit = models.CharField(
        max_length=50,
        null=True,
        blank=True,
        help_text="Capacity unit, e.g. ul, ml, L, mm, um, cm2, gal.",
    )

    description = models.TextField(
        null=True,
        blank=True,
        help_text="Optional free-text description for extra metadata not captured elsewhere.",
    )

    default_cost = models.DecimalField(
        max_digits=12,
        decimal_places=2,
        null=True,
        blank=True,
        help_text="Optional default/reference cost for one order unit or stock unit.",
    )

    lifetime_days = models.PositiveIntegerField(
        null=True,
        blank=True,
        help_text="Optional default lifetime in days, mainly useful for equipment or reusable items.",
    )

    serial_number = models.CharField(
        max_length=255,
        null=True,
        blank=True,
        help_text="Optional serial number. Usually relevant only for devices/accessories, not consumables.",
    )

    order_number = models.CharField(
        max_length=255,
        null=True,
        blank=True,
        help_text="Optional internal order reference.",
    )

    is_active = models.BooleanField(
        default=True,
        help_text="Soft flag to hide discontinued materials without deleting history.",
    )

    class Meta:
        ordering = ["product_name"]
        verbose_name = "Material master"
        verbose_name_plural = "Material master records"

    def __str__(self):
        return self.product_name


class MaterialUnit(models.Model):
    """
    Defines which counting / packaging units are valid for a material
    and how they convert to the base unit.

    This lets one material support both scenarios:
    - counted in boxes
    - counted in individual pieces

    Example: 200ul MCA96 Tips
    - tip       -> base unit, factor = 1
    - tip box   -> stock unit, factor = 96
    - carton    -> order unit, factor = 960

    Example: Gloves M
    - glove     -> base unit, factor = 1
    - box       -> stock unit, factor = 100

    Example: Prime spare power cable
    - piece     -> base unit, stock unit, order unit, factor = 1
    """

    material = models.ForeignKey(
        MaterialMaster,
        on_delete=models.CASCADE,
        related_name="units",
        help_text="Material to which this unit definition belongs.",
    )

    unit = models.ForeignKey(
        UnitOfMeasure,
        on_delete=models.PROTECT, # PROTECT means models.PROTECT prevents deletion of the referenced object if it is still being used by another model.
        related_name="material_units",
        help_text="Unit name, e.g. piece, box, bag, rack, plate, bottle.",
    )

    is_base_unit = models.BooleanField(
        default=False,
        help_text="Marks the smallest logical unit for this material, e.g. tip, glove, tube, piece.",
    )

    is_stock_unit = models.BooleanField(
        default=False,
        help_text="Marks the unit usually used in inventory stock records, e.g. box, bag, plate, piece.",
    )

    is_order_unit = models.BooleanField(
        default=False,
        help_text="Marks the unit usually used when ordering from suppliers.",
    )
    # NOTE:
    # This field stores how many base units are contained inside one packaging unit.
    #
    # Example conversions:
    # - 1 tip       = 1 base unit
    # - 1 tip box   = 96 tips
    # - 1 carton    = 960 tips
    # - 1 glove box = 100 gloves
    #
    # We intentionally use DecimalField instead of IntegerField because:
    # - inventory quantities can be fractional (e.g. 0.8 open boxes remaining)
    # - Excel source data already contains fractional packaging values
    # - some packaging conversions are not whole numbers
    # - float values would introduce rounding errors (bad for inventory/accounting)
    #
    # Example:
    # 2.041667 boxes × 96 tips = 196 tips exactly
    #
    # decimal_places=6 matches the precision observed in the original Excel inventory
    # and ensures conversions remain accurate during stock tracking and billing.

    base_units_per_unit = models.DecimalField(
        max_digits=12,
        decimal_places=6,
        default=1,
        help_text=(
            "How many base units are contained in one of this unit. "
            "Examples: tip=1, tip box=96, glove box=100, carton=960."
        ),
    )

    notes = models.TextField(
        null=True,
        blank=True,
        help_text=(
            "Optional explanation for unusual packaging from our original excel, "
            "e.g. '24 racks in 2er blister', 'single reservoirs', '1 open box'."
        ),
    )

    class Meta:
        unique_together = ("material", "unit")
        ordering = ["material__product_name", "unit__name"]
        verbose_name = "Material unit"
        verbose_name_plural = "Material units"

    def __str__(self):
        return f"{self.material.product_name} [{self.unit.name}]"