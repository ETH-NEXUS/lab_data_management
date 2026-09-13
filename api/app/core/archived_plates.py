"""
Keeps archived library plates from being changed through the API.

The lab archives library plates that were physically removed from storage. Their data
stays visible, but the plates and their wells may not change any more. Everything else
keeps working: experiments that used such a plate can still be read and measured, and a
plate that is not a library plate is never blocked. The admin works on the models
directly and is not affected, so a mistake can still be corrected there.
"""

from rest_framework.exceptions import ValidationError


def ensure_plate_can_be_changed(plate):
    """
    Refuse the request with a 400 response when the plate is an archived library plate.
    Nothing happens for `None`, so callers can pass an optional plate as it is.
    Response data example:
    {"detail": "Plate Drug03_I is archived and can no longer be changed."}
    """
    if plate is not None and plate.archived and plate.library_id is not None:
        raise ValidationError(
            {
                "detail": f"Plate {plate.barcode} is archived and can no longer be changed."
            }
        )
