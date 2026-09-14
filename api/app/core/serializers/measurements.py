"""
Measurements of wells and their features.
"""

from rest_framework import serializers

from ..models import Measurement, MeasurementFeature


class MeasurementFeatureSerializer(serializers.ModelSerializer):
    class Meta:
        model = MeasurementFeature
        fields = "__all__"


class MeasurementSerializer(serializers.ModelSerializer):
    feature = MeasurementFeatureSerializer()

    class Meta:
        model = Measurement
        fields = "__all__"

    def to_representation(self, instance):
        representation = super().to_representation(instance)
        representation["measured_at"] = instance.measured_at.isoformat().split("+")[0]

        return representation
