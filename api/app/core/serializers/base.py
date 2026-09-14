"""
Base serializers shared by other serializers of the core app.
"""

from rest_framework import serializers


class UndefinedAffineModelSerializer(serializers.ModelSerializer):
    """
    A serializer that takes care about undefined values in the request
    data and converts them to None.
    """

    def to_internal_value(self, data):
        for key, value in data.items():
            if value == "undefined":
                data[key] = None
        return super().to_internal_value(data)
