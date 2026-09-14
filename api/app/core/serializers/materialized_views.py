"""
Read-only details of wells, plates and experiments from the materialized views.
"""

from rest_framework import serializers

from ..models import ExperimentDetail, PlateDetail, WellDetail


class WellDetailSerializer(serializers.ModelSerializer):
    class Meta:
        model = WellDetail
        fields = "__all__"


class PlateDetailSerializer(serializers.ModelSerializer):
    empty = {}

    class Meta:
        model = PlateDetail
        fields = "__all__"


class ExperimentDetailSerializer(serializers.ModelSerializer):
    empty = {}

    class Meta:
        model = ExperimentDetail
        fields = "__all__"
