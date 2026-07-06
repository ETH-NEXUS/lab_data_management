from django.contrib.auth.models import User
from django.db import models

from core.models import Experiment, Project

from .static_models import MaterialUnit


class InventoryChangeRecord(models.Model):
    """
    Immutable audit row for one inventory-related action.

    This model stores a generic history stream for stock, order, and usage actions.

    Example rows:
    - {
        "performed_action": "stock_created",
        "performed_by": 3,
        "quantity_delta": "2.000000",
        "notes": "Created from stock create endpoint"
      }
    - {
        "performed_action": "stock_archived",
        "performed_by": 7,
        "notes": "Archived from inventory stock detail view"
      }
    """

    ACTION_STOCK_CREATED = "stock_created"
    ACTION_STOCK_UPDATED = "stock_updated"
    ACTION_STOCK_DELETED = "stock_deleted"
    ACTION_STOCK_ARCHIVED = "stock_archived"
    ACTION_STOCK_RESTORED = "stock_restored"
    ACTION_STOCK_FAVORITED = "stock_favorited"
    ACTION_STOCK_UNFAVORITED = "stock_unfavorited"
    ACTION_ORDER_CREATED = "order_created"
    ACTION_ORDER_UPDATED = "order_updated"
    ACTION_ORDER_DELETED = "order_deleted"
    ACTION_USAGE_CREATED = "usage_created"
    ACTION_USAGE_UPDATED = "usage_updated"
    ACTION_USAGE_DELETED = "usage_deleted"

    ACTION_CHOICES = [
        (ACTION_STOCK_CREATED, "Stock created"),
        (ACTION_STOCK_UPDATED, "Stock updated"),
        (ACTION_STOCK_DELETED, "Stock deleted"),
        (ACTION_STOCK_ARCHIVED, "Stock archived"),
        (ACTION_STOCK_RESTORED, "Stock restored"),
        (ACTION_STOCK_FAVORITED, "Stock favorited"),
        (ACTION_STOCK_UNFAVORITED, "Stock unfavorited"),
        (ACTION_ORDER_CREATED, "Order created"),
        (ACTION_ORDER_UPDATED, "Order updated"),
        (ACTION_ORDER_DELETED, "Order deleted"),
        (ACTION_USAGE_CREATED, "Usage created"),
        (ACTION_USAGE_UPDATED, "Usage updated"),
        (ACTION_USAGE_DELETED, "Usage deleted"),
    ]

    performed_action = models.CharField(
        max_length=100,
        choices=ACTION_CHOICES,
        help_text="Type of inventory action that was performed.",
    )

    performed_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="inventory_change_records",
        help_text="Authenticated user who performed the action if available.",
    )

    performed_at = models.DateTimeField(
        auto_now_add=True,
        help_text="When the action was recorded.",
    )

    inventory_stock = models.ForeignKey(
        "InventoryStock",
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="change_records",
        help_text="Related stock row if the action affected stock.",
    )

    order = models.ForeignKey(
        "Order",
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="change_records",
        help_text="Related order row if the action affected an order.",
    )

    material_usage = models.ForeignKey(
        "MaterialUsage",
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="change_records",
        help_text="Related usage row if the action affected a material usage entry.",
    )

    project = models.ForeignKey(
        Project,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="inventory_change_records",
        help_text="Optional related project for the recorded action.",
    )

    experiment = models.ForeignKey(
        Experiment,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="inventory_change_records",
        help_text="Optional related experiment for the recorded action.",
    )

    quantity_delta = models.DecimalField(
        max_digits=12,
        decimal_places=6,
        null=True,
        blank=True,
        help_text=(
            "Optional quantity change for the action. "
            "Example: {\"quantity_delta\": \"-1.000000\"} for a checkout."
        ),
    )

    quantity_unit = models.ForeignKey(
        MaterialUnit,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="inventory_change_records",
        help_text="Optional unit used with quantity_delta.",
    )

    notes = models.TextField(
        null=True,
        blank=True,
        help_text="Optional free-text explanation for the action.",
    )

    class Meta:
        ordering = ["-performed_at", "-id"]
        verbose_name = "Inventory change record"
        verbose_name_plural = "Inventory change records"

    def __str__(self):
        return f"{self.performed_action} @ {self.performed_at}"
