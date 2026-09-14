"""
Plates, their dimensions and plate templates.
"""

from rest_framework import serializers

from compoundlib.models import CompoundLibrary
from platetemplate.models import PlateTemplate
from ..models import Plate, PlateDetail, PlateDimension, WellDetail
from .materialized_views import PlateDetailSerializer, WellDetailSerializer


class SimplePlateSerializer(serializers.ModelSerializer):
    dimension = serializers.SlugRelatedField(read_only=True, slug_field="name")
    library = serializers.SlugRelatedField(read_only=True, slug_field="name")
    measurement_labels = serializers.SerializerMethodField()

    def get_measurement_labels(self, plate: Plate):
        try:
            plate_details = PlateDetail.objects.get(pk=plate.id)
            return plate_details.measurement_labels
        except PlateDetail.DoesNotExist:
            return PlateDetailSerializer.empty

    class Meta:
        model = Plate
        fields = (
            "id",
            "barcode",
            "dimension",
            "library",
            "measurement_labels",
            "archived",
            "status",
        )


class PlateDimensionSerializer(serializers.ModelSerializer):
    class Meta:
        model = PlateDimension
        fields = "__all__"
        extra_kwargs = {
            "id": {
                "read_only": False,
                "required": False,
            },
        }


class PlateSerializer(serializers.ModelSerializer):
    dimension = PlateDimensionSerializer(required=False, allow_null=True)
    details = serializers.SerializerMethodField()
    wells = serializers.SerializerMethodField()

    def get_details(self, plate: Plate):
        try:
            plate_details = PlateDetail.objects.get(pk=plate.id)
            return PlateDetailSerializer(plate_details).data
        except PlateDetail.DoesNotExist:
            return PlateDetailSerializer.empty

    def get_wells(self, plate: Plate):
        wells = WellDetail.objects.filter(plate_id=plate.id).order_by("position")
        return WellDetailSerializer(wells, many=True).data

    def update(self, plate: Plate, validated_data):
        if "dimension" in validated_data:
            dimension = validated_data.pop("dimension")
            dimension_id = dimension.get("id")
            if dimension_id:
                plate.dimension = PlateDimension.objects.get(pk=dimension_id)
        if "library" in validated_data:
            library = validated_data.pop("library")
            library_id = library.get("id")
            if library_id:
                plate.library = CompoundLibrary.objects.get(pk=library_id)
        plate.save()

        return plate

    class Meta:
        model = Plate
        fields = "__all__"


class SimplePlateTemplateSerializer(serializers.ModelSerializer):
    category = serializers.SlugRelatedField(slug_field="name", read_only=True)

    class Meta:
        model = PlateTemplate
        fields = ("id", "name", "category")
