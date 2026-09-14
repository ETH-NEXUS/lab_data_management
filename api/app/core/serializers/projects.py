"""
Projects, with their experiments and plates.
"""

import random
import string

from rest_framework import serializers

from ..models import Project
from .experiments import ExperimentSerializer
from .plates import PlateSerializer


class ProjectSerializer(serializers.ModelSerializer):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if self.instance:
            if isinstance(self.instance, list):
                for instance in self.instance:
                    if instance.harvest_id is not None:
                        self.fields["name"].read_only = True
                        break
            else:
                if self.instance.harvest_id is not None:
                    self.fields["name"].read_only = True

    experiments = ExperimentSerializer(many=True, required=False, allow_null=True)
    plates = PlateSerializer(many=True, required=False, allow_null=True)

    class Meta:
        model = Project
        fields = "__all__"

    def to_representation(self, instance):
        representation = super().to_representation(instance)

        request = self.context.get("request")
        if request and request.user.username == "demo":
            random_suffix = "".join(
                random.choices(string.ascii_letters + string.digits, k=10)
            )
            representation["name"] = f"Project_{random_suffix}"

        return representation


class SimpleProjectSerializer(serializers.ModelSerializer):
    class Meta:
        model = Project
        fields = ("id", "name")
