"""
Starting the background correction of one plate from the plate page.
"""

from django.shortcuts import get_object_or_404
from rest_framework import serializers
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response

from background_correction.calculation import METHODS
from background_correction.correction import correct_plate
from core.models import Plate


class CorrectionSettingsSerializer(serializers.Serializer):
    label = serializers.CharField()
    reference_type = serializers.CharField()
    method = serializers.ChoiceField(choices=list(METHODS))


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def correct_plate_background(request, plate_id: int):
    """
    Saves the corrected measurement of the plate and returns its label.

    Accepted data example:
    {"label": "Lum_CTG", "reference_type": "N1", "method": "median"}
    Returned data example:
    {"label": "Lum_CTG_bc_N1_median"}
    """
    plate = get_object_or_404(Plate, id=plate_id)
    settings = CorrectionSettingsSerializer(data=request.data)
    settings.is_valid(raise_exception=True)

    new_label = correct_plate(
        plate,
        settings.validated_data["label"],
        settings.validated_data["reference_type"],
        settings.validated_data["method"],
    )
    return Response({"label": new_label})
