from core.models import Project
from inventory.dynamic_models import Order


def get_or_create_named(model_class, name):
    """
    Generic helper for lookup tables with a 'name' field.
    """
    if not name:
        return None

    obj, _ = model_class.objects.get_or_create(name=name)
    return obj


def normalize_order_status(value):
    """
    Normalize free-text excel order status to model choices.
    """
    value = (value or "").strip().lower()

    if "tentative" in value:
        return Order.STATUS_TENTATIVE

    if "arrived" in value:
        return Order.STATUS_PRODUCT_ARRIVED

    return Order.STATUS_ORDERED


def find_project_by_name(project_name):
    """
    Try to match a project by exact name.
    """
    if not project_name:
        return None

    return Project.objects.filter(name=project_name).first()