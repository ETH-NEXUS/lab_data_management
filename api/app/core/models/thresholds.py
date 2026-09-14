"""
The threshold below which a library well counts as running low.
"""

from django.core.validators import MaxValueValidator, MinValueValidator
from django.db import models


class Threshold(models.Model):
    """
    The limits below which a library well counts as running low.
    There is exactly one of these for the whole application; always read it
    through `Threshold.current()`.
    """

    dmso = models.FloatField(
        default=80,
        validators=[MinValueValidator(0), MaxValueValidator(100)],
        help_text="In percent. A well with less DMSO than this is marked.",
    )
    amount = models.FloatField(
        default=2.5,
        validators=[MinValueValidator(0)],
        help_text="In microliter (µL), the volume left in the source well. "
        "A well with less than this is marked.",
    )

    @classmethod
    def current(cls) -> "Threshold":
        """
        Return the one threshold, creating it with the defaults if it is missing.
        Returned data example:
        {"id": 1, "amount": 2.5, "dmso": 80}
        """
        threshold = cls.objects.order_by("id").first()
        if threshold is None:
            # A fixed primary key makes two requests that both find the table
            # empty collide instead of creating two rows; get_or_create then
            # simply reads the row the other request created.
            threshold, _ = cls.objects.get_or_create(pk=1)
        return threshold
