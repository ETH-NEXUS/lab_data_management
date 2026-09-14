"""
Experiments and their barcode specifications.
"""

from collections import defaultdict

from rest_framework import serializers

from ..models import BarcodeSpecification, Experiment, ExperimentDetail, PlateDetail
from .materialized_views import ExperimentDetailSerializer
from .plates import SimplePlateSerializer


class BarcodeSpecificationSerializer(serializers.ModelSerializer):
    class Meta:
        model = BarcodeSpecification
        fields = "__all__"


class ExperimentSerializer(serializers.ModelSerializer):
    plates = SimplePlateSerializer(many=True, required=False, allow_null=True)
    barcode_specifications = BarcodeSpecificationSerializer(
        many=True, required=False, allow_null=True
    )
    available_measurement_labels = serializers.SerializerMethodField()

    details = serializers.SerializerMethodField()

    def get_details(self, experiment: Experiment):
        try:
            experiment_details = ExperimentDetail.objects.get(pk=experiment.id)
            return ExperimentDetailSerializer(experiment_details).data
        except ExperimentDetail.DoesNotExist:
            return ExperimentDetailSerializer.empty

    # returns only those labels that are available for all plates in the experiment
    def get_available_measurement_labels(self, experiment: Experiment):
        labels_count = defaultdict(int)
        labels = []
        for plate in experiment.plates.all():
            try:
                plate_details = PlateDetail.objects.get(pk=plate.id)
                labels.extend(plate_details.measurement_labels)
            except PlateDetail.DoesNotExist:
                pass
        for label in labels:
            labels_count[label] += 1
        return [
            label
            for label, count in labels_count.items()
            if count == len(experiment.plates.all())
        ]

    class Meta:
        model = Experiment
        fields = "__all__"


class SimpleExperimentSerializer(serializers.ModelSerializer):
    class Meta:
        model = Experiment
        fields = ("id", "name")
