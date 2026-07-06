from .history_models import InventoryChangeRecord


def record_inventory_action(
    *,
    performed_action,
    performed_by=None,
    inventory_stock=None,
    order=None,
    material_usage=None,
    project=None,
    experiment=None,
    quantity_delta=None,
    quantity_unit=None,
    notes=None,
):
    """
    Stores one immutable inventory history row.

    Accepted example:
    - {
        "performed_action": "stock_created",
        "performed_by": <User>,
        "inventory_stock": <InventoryStock>,
        "quantity_delta": "2.000000",
        "quantity_unit": <MaterialUnit>,
        "notes": "Created from stock endpoint"
      }

    Returned example:
    - {
        "id": 14,
        "performed_action": "stock_created"
      }
    """

    change_record = InventoryChangeRecord.objects.create(
        performed_action=performed_action,
        performed_by=performed_by,
        inventory_stock=inventory_stock,
        order=order,
        material_usage=material_usage,
        project=project,
        experiment=experiment,
        quantity_delta=quantity_delta,
        quantity_unit=quantity_unit,
        notes=notes,
    )

    return change_record
