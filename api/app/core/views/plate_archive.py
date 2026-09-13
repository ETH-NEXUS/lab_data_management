"""
Archiving and unarchiving a plate, for the button on the plate page.
"""

from rest_framework import status
from rest_framework.decorators import action
from rest_framework.generics import get_object_or_404
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response

from ..models import Plate


class PlateArchiveMixin:
    """
    Adds POST /api/plates/<id>/archive/ to the plate viewset.

    This is a separate action on purpose: the general plate update only saves the
    dimension and the library, and changing that shared behaviour would affect every
    other way a plate is edited.
    """

    # POST, so DRF checks the CSRF token. Every logged in user may archive a plate.
    @action(detail=True, methods=["post"], permission_classes=[IsAuthenticated])
    def archive(self, request, pk=None):
        """
        Archives or unarchives one plate.
        Accepted data example:
        {"archived": true}
        Returned data example:
        {"id": 7, "barcode": "Drug01_E", "archived": true}
        """
        archived = request.data.get("archived")
        if not isinstance(archived, bool):
            return Response(
                {
                    "archived": [
                        "Send true to archive the plate or false to unarchive it."
                    ]
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        # Looked up directly: the queryset of the plate viewset loads every well of the
        # plate with its withdrawals and measurements, which this action does not need.
        plate = get_object_or_404(Plate, pk=pk)
        plate.archived = archived
        plate.save(update_fields=["archived", "modified_at"])

        return Response(
            {"id": plate.id, "barcode": plate.barcode, "archived": plate.archived}
        )
