from django.db import models
from core.models import Project, Experiment
from .static_models import MaterialMaster, MaterialUnit
from django.contrib.auth.models import User

"""
These dymanic models show:

	•	where it is
	•	how much is there
	•	what batch it belongs to
	•	whether it is low stock
	•	who ordered it
	•	how much was used by a project/experiment

"""

class Room(models.Model):
    """
    Physical room / storage area.

    Examples from Excel:
    - C75
    - C41
    - C21
    - C51
    """
    name = models.CharField(
        max_length=255,
        unique=True,
        help_text="Room name, e.g. C75, C41, C21.",
    )

    class Meta:
        ordering = ["name"]
        verbose_name = "Room"
        verbose_name_plural = "Rooms"

    def __str__(self):
        return self.name


class Sector(models.Model):
    """
    Shelf / cupboard / freezer sector inside a room.

    Examples from the lab's excel:
    - 11.2
    - 7.1
    - 4.3
    - freezer
    - 3.1;3.2
    """
    name = models.CharField(
        max_length=255,
        help_text="Sector / shelf label inside the room.",
    )
    room = models.ForeignKey(
        Room,
        on_delete=models.CASCADE,
        related_name="sectors",
    )

    class Meta:
        ordering = ["room__name", "name"]
        unique_together = ("room", "name")
        verbose_name = "Sector"
        verbose_name_plural = "Sectors"

    def __str__(self):
        return f"{self.room} / {self.name}"


class InventoryStock(models.Model):
    """
    One row = one current stock record for a material in a location.
    This is the dynamic equivalent of one excel inventory row.

    It stores:
    - which material
    - where it is
    - lot / expiry for this batch
    - how much is currently available
    - in which unit it is counted
    - threshold for low stock

    Examples:
    - 1.2 boxes of a PCR plate lot in C75 / 5.4
    - 3 boxes of gloves in C75 / 6.3
    - 1 piece of a spare cable in C41 / 1.3
    """

    material = models.ForeignKey(
        MaterialMaster,
        on_delete=models.CASCADE,
        related_name="stock_entries",
        help_text="Material master record for this stock item.",
    )

    sector = models.ForeignKey(
        Sector,
        on_delete=models.CASCADE,
        related_name="stock_entries",
        help_text="Physical storage location of this stock entry.",
    )

    stock_unit = models.ForeignKey(
        MaterialUnit,
        on_delete=models.PROTECT,
        related_name="stock_entries",
        help_text="Unit in which this stock entry is counted, e.g. box, piece, bag.",
    )

    quantity = models.DecimalField(
        max_digits=12,
        decimal_places=6,
        default=0,
        help_text=(
            "Current stock quantity in the selected stock unit. "
            "Can be fractional for open boxes, e.g. 0.8, 1.25, 2.041667."
        ),
    )

    # Minimum desired stock threshold expressed in the same unit as `quantity`.
    #
    # This is intentionally a DecimalField (not IntegerField) because:
    # - stock quantities may be fractional (e.g. 0.8 open boxes remaining)
    # - lab can define reorder thresholds below one full package (e.g. 0.5 box)
    # - comparisons with `quantity` must remain precise
    # - FloatField would introduce rounding errors in inventory logic
    #
    # Example:
    # quantity = 0.8 boxes
    # minimum_quantity = 1 box
    # → inventory_status = "low"

    minimum_quantity = models.DecimalField(
        max_digits=12,
        decimal_places=6,
        default=0,
        help_text=(
            "Minimum desired quantity in the same unit as 'quantity'. "
            "Used to derive low-stock status."
        ),
    )

    lot_number = models.CharField(
        max_length=255,
        null=True,
        blank=True,
        help_text="Lot / batch number for this specific stock entry.",
    )

    expiry_date = models.DateField(
        null=True,
        blank=True,
        help_text="Expiry date for this specific lot / batch if applicable.",
    )

    is_favorite = models.BooleanField(
        default=False,
        help_text="Optional flag if we want to mark favorite / pinned stock items.",
    )
    is_archived = models.BooleanField(
        default=False,
        help_text="Optional flag if we want to hide archived stock items.",
    )

    notes = models.TextField(
        null=True,
        blank=True,
        help_text=(
            "Free-text stock notes, e.g. '1 open box', 'outside box', "
            "'all open', 'very old', 'single reservoirs'."
        ),
    )

    created_at = models.DateTimeField(
        auto_now_add=True,
        help_text="When this stock entry was created in the system.",
    )

    updated_at = models.DateTimeField(
        auto_now=True,
        help_text="When this stock entry was last updated.",
    )

    class Meta:
        ordering = ["material__product_name", "sector__room__name", "sector__name"]
        verbose_name = "Inventory stock"
        verbose_name_plural = "Inventory stock"

    def __str__(self):
        return f"{self.material} @ {self.sector}"

    @property
    def quantity_in_base_units(self):
        """
        Converts the current quantity into the material's base unit.
        Example:
        - 2 boxes * 96 tips per box = 192 tips
        """
        return self.quantity * self.stock_unit.base_units_per_unit

    @property
    def minimum_quantity_in_base_units(self):
        """
        Converts the minimum quantity into the material's base unit.
        """
        return self.minimum_quantity * self.stock_unit.base_units_per_unit

    @property
    def inventory_status(self):
        """
        Equivalent of the Excel computed status:
        - low if current quantity < minimum quantity
        - otherwise in stock
        """
        if self.quantity < self.minimum_quantity:
            return "low"
        return "in stock"


class Order(models.Model):
    """
    Purchase / ordering record for a material.

    This is not stock itself.
    It records an ordering event and can be linked to a project.

    Examples of statuses from Excel:
    - ordered
    - tentative
    - product arrived
    """

    STATUS_ORDERED = "ordered"
    STATUS_TENTATIVE = "tentative"
    STATUS_PRODUCT_ARRIVED = "product_arrived"

    STATUS_CHOICES = [
        (STATUS_ORDERED, "Ordered"),
        (STATUS_TENTATIVE, "Tentative"),
        (STATUS_PRODUCT_ARRIVED, "Product arrived"),
    ]

    material = models.ForeignKey(
        MaterialMaster,
        on_delete=models.CASCADE,
        related_name="orders",
        help_text="Material being ordered.",
    )

    order_unit = models.ForeignKey(
        MaterialUnit,
        on_delete=models.PROTECT,
        related_name="orders",
        help_text="Unit in which the order was placed, e.g. box, carton, piece.",
    )

    amount = models.DecimalField(
        max_digits=12,
        decimal_places=6,
        help_text="Ordered quantity in the selected order unit.",
    )

    order_date = models.DateTimeField(
        help_text="Date/time when the order was placed or recorded.",
    )

    status = models.CharField(
        max_length=50,
        choices=STATUS_CHOICES,
        default=STATUS_ORDERED,
        help_text="Current order status.",
    )

    who_ordered = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="material_orders",
        help_text="User who placed the order.", )

    project = models.ForeignKey(
        Project,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="material_orders",
        help_text="Optional project to which this order belongs.",
    )

    notes = models.TextField(
        null=True,
        blank=True,
        help_text="Free-text notes about the order.",
    )

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        ordering = ["-order_date"]
        verbose_name = "Order"
        verbose_name_plural = "Orders"

    def __str__(self):
        return f"Order: {self.material}"


class MaterialUsage(models.Model):
    """
    Material consumption / assignment to a project or experiment.


    usage is stored with an explicit unit, so the lab can record:
    - 2 boxes
    - 12 plates
    - 120 tips
    - 1 piece

    This is what enables charging projects for individual pieces
    even when stock is stored in boxes.
    """

    project = models.ForeignKey(
        Project,
        on_delete=models.CASCADE,
        related_name="material_usages",
        help_text="Project to which the material usage belongs.",
    )

    experiment = models.ForeignKey(
        Experiment,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="material_usages",
        help_text="Optional experiment within the project.",
    )

    inventory_stock = models.ForeignKey(
        InventoryStock,
        on_delete=models.PROTECT,
        related_name="usages",
        help_text="Stock entry from which the material was consumed.",
    )

    usage_unit = models.ForeignKey(
        MaterialUnit,
        on_delete=models.PROTECT,
        related_name="usage_entries",
        help_text="Unit used for this consumption entry, e.g. piece, tip, box, plate.",
    )

    quantity_used = models.DecimalField(
        max_digits=12,
        decimal_places=6,
        help_text="Consumed quantity in the selected usage unit.",
    )

    used_at = models.DateTimeField(
        auto_now_add=True,
        help_text="When this material usage was recorded.",
    )

    notes = models.TextField(
        null=True,
        blank=True,
        help_text="Optional free-text notes about the usage.",
    )

    class Meta:
        ordering = ["-used_at"]
        verbose_name = "Material usage"
        verbose_name_plural = "Material usages"

    def __str__(self):
        return f"{self.project} - {self.inventory_stock.material}"

    @property
    def quantity_used_in_base_units(self):
        """
        Converts the used quantity into the material's base unit.

        Example:
        - 2 boxes * 96 tips per box = 192 tips
        """
        return self.quantity_used * self.usage_unit.base_units_per_unit
