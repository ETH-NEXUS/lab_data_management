from django.utils import timezone


def get_current_date_safe():
    """
    Returns today's date safely for both timezone-aware and naive datetime setups.

    Returned data examples:
    - aware setup: `timezone.localdate(now)`
    - naive setup: `timezone.now().date()`
    """
    now_value = timezone.now()

    if timezone.is_naive(now_value):
        return now_value.date()

    return timezone.localdate(now_value)
