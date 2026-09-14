"""
The threshold below which a library well counts as running low.
"""

from rest_framework import serializers

from ..models import Threshold


class ThresholdSerializer(serializers.ModelSerializer):
    class Meta:
        model = Threshold
        fields = "__all__"
