"""
Starting the calculations of one plate from the plate page: background
correction, log10 and %Activity.
"""

from django.shortcuts import get_object_or_404
from rest_framework import serializers
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response

from core.models import Plate
from plate_calculations.correction import METHODS, correct_plate
from plate_calculations.log_transform import log10_of_plate
from plate_calculations.percent_activity import activity_of_plate


class CorrectionSettingsSerializer(serializers.Serializer):
    # Labels and well types are looked up exactly: some labels end with a space
    label = serializers.CharField(trim_whitespace=False)
    reference_type = serializers.CharField(trim_whitespace=False)
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


class Log10SettingsSerializer(serializers.Serializer):
    # Labels are looked up exactly: some labels end with a space
    label = serializers.CharField(trim_whitespace=False)


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def log10_plate_measurement(request, plate_id: int):
    """
    Saves the log10 of a measurement of the plate; wells with a value of 0 or
    below are left empty and counted.

    Accepted data example:
    {"label": "Lum1"}
    Returned data example:
    {"label": "Lum1_log10", "skipped": 2}
    """
    plate = get_object_or_404(Plate, id=plate_id)
    settings = Log10SettingsSerializer(data=request.data)
    settings.is_valid(raise_exception=True)

    new_label, skipped = log10_of_plate(plate, settings.validated_data["label"])
    return Response({"label": new_label, "skipped": skipped})


class ActivitySettingsSerializer(serializers.Serializer):
    # Labels and well types are looked up exactly: some labels end with a space
    label = serializers.CharField(trim_whitespace=False)
    negative_type = serializers.CharField(trim_whitespace=False)
    positive_type = serializers.CharField(trim_whitespace=False)

    def validate(self, data):
        if data["negative_type"] == data["positive_type"]:
            raise serializers.ValidationError(
                "The negative and the positive control must be different well types."
            )
        return data


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def percent_activity_of_plate(request, plate_id: int):
    """
    Saves the %Activity of a measurement of the plate between its controls.

    Accepted data example:
    {"label": "Lum1", "negative_type": "N", "positive_type": "P"}
    Returned data example:
    {"label": "Lum1_activity_N_P"}
    """
    plate = get_object_or_404(Plate, id=plate_id)
    settings = ActivitySettingsSerializer(data=request.data)
    settings.is_valid(raise_exception=True)

    new_label = activity_of_plate(
        plate,
        settings.validated_data["label"],
        settings.validated_data["negative_type"],
        settings.validated_data["positive_type"],
    )
    return Response({"label": new_label})
