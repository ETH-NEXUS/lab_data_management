"""
Starting the calculations of one plate from the plate page: background
correction, log10 and the normalization (%Inhibition, %Activity).
"""

from django.shortcuts import get_object_or_404
from rest_framework import serializers
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response

from core.models import Plate
from core.utils.plates.archive_guard import ensure_plate_can_be_changed
from plate_calculations.correction import METHODS, correct_plate
from plate_calculations.log_transform import log10_of_plate
from plate_calculations.normalization import normalize_plate


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
    ensure_plate_can_be_changed(plate)
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
    ensure_plate_can_be_changed(plate)
    settings = Log10SettingsSerializer(data=request.data)
    settings.is_valid(raise_exception=True)

    new_label, skipped = log10_of_plate(plate, settings.validated_data["label"])
    return Response({"label": new_label, "skipped": skipped})


class NormalizationSettingsSerializer(serializers.Serializer):
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
def normalize_plate_measurement(request, plate_id: int):
    """
    Saves the %Inhibition and %Activity of a measurement of the plate between its
    controls; wells with a value of 0 or below are left empty and counted.

    Accepted data example:
    {"label": "Lum1", "negative_type": "N", "positive_type": "P"}
    Returned data example (label: the %Inhibition, shown first):
    {"label": "Lum1_inhibition_N_P", "activity_label": "Lum1_activity_N_P", "skipped": 0}
    """
    plate = get_object_or_404(Plate, id=plate_id)
    ensure_plate_can_be_changed(plate)
    settings = NormalizationSettingsSerializer(data=request.data)
    settings.is_valid(raise_exception=True)

    inhibition_label, activity_label, skipped = normalize_plate(
        plate,
        settings.validated_data["label"],
        settings.validated_data["negative_type"],
        settings.validated_data["positive_type"],
    )
    return Response(
        {
            "label": inhibition_label,
            "activity_label": activity_label,
            "skipped": skipped,
        }
    )
