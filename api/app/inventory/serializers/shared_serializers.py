from django.contrib.auth import get_user_model
from rest_framework import serializers

User = get_user_model()


class SimpleNameSerializer(serializers.Serializer):
    """
    Small reusable serializer for objects with:
    - id
    - name
    - string label
    """
    id = serializers.IntegerField(read_only=True)
    name = serializers.CharField(read_only=True)
    label = serializers.SerializerMethodField()

    def get_label(self, obj):
        return str(obj)


class SimpleObjectLabelSerializer(serializers.Serializer):
    """
    Small reusable serializer for objects where we only want:
    - id
    - string label
    """
    id = serializers.IntegerField(read_only=True)
    label = serializers.SerializerMethodField()

    def get_label(self, obj):
        return str(obj)


class UserSummarySerializer(serializers.Serializer):
    """
    Small user serializer for UI display.
    """
    id = serializers.IntegerField(read_only=True)
    username = serializers.CharField(read_only=True)
    full_name = serializers.SerializerMethodField()
    label = serializers.SerializerMethodField()

    def get_full_name(self, obj):
        full_name = obj.get_full_name()
        return full_name if full_name else obj.username

    def get_label(self, obj):
        full_name = obj.get_full_name()
        return full_name if full_name else obj.username