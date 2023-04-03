from compoundlib.models import CompoundLibrary, Compound
from rest_framework import serializers

from .models import (
    Well,
    Plate,
    Measurement,
    MeasurementFeature,
    PlateDimension,
    Project,
    Experiment,
    WellCompound,
    WellWithdrawal,
    PlateMapping,
    PlateTemplate,
    MappingError,
    WellType,
    BarcodeSpecification,
    MappingList,
    PlateDetail,
    WellDetail,
)


class UndefinedAffineModelSerializer(serializers.ModelSerializer):
    """
    A serializer that takes care about undefined values in the request
    data and converts them to None.
    """

    def to_internal_value(self, data):
        for key, value in data.items():
            print("value", value)
            if value == "undefined":
                data[key] = None
        return super().to_internal_value(data)


class CompoundSerializer(serializers.ModelSerializer):
    class Meta:
        model = Compound
        fields = ("name", "identifier", "structure")


class WellCompoundSerializer(serializers.ModelSerializer):
    name = serializers.SlugRelatedField(
        slug_field="name", source="compound", read_only=True
    )
    identifier = serializers.SlugRelatedField(
        slug_field="identifier", source="compound", read_only=True
    )
    structure = serializers.SlugRelatedField(
        slug_field="structure", source="compound", read_only=True
    )

    class Meta:
        model = WellCompound
        fields = "__all__"


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


class SimplePlateSerializer(serializers.ModelSerializer):
    dimension = serializers.SlugRelatedField(read_only=True, slug_field="name")
    library = serializers.SlugRelatedField(read_only=True, slug_field="name")

    class Meta:
        model = Plate
        fields = ("id", "barcode", "dimension", "library")


class WellPlateSerializer(serializers.ModelSerializer):
    plate = SimplePlateSerializer()
    hr_position = serializers.ReadOnlyField()
    amount = serializers.ReadOnlyField()
    mixture = serializers.SerializerMethodField()

    def get_mixture(self, well: Well):
        """If this is a well with multiple compounds"""
        return well.compounds.count() > 1

    class Meta:
        model = Well
        exclude = ("compounds",)


class WellWithdrawalSerializer(serializers.ModelSerializer):
    well = WellPlateSerializer()
    target_well = WellPlateSerializer()

    class Meta:
        model = WellWithdrawal
        fields = "__all__"


class WellTypeSerializer(serializers.ModelSerializer):
    class Meta:
        model = WellType
        fields = "__all__"


class WellDetailSerializer(serializers.ModelSerializer):
    class Meta:
        model = WellDetail
        fields = "__all__"


class WellSerializer(serializers.ModelSerializer):
    hr_position = serializers.ReadOnlyField()
    compounds = WellCompoundSerializer(
        many=True, required=False, allow_null=True, source="well_compounds"
    )
    withdrawals = WellWithdrawalSerializer(many=True, required=False, allow_null=True)
    donors = WellWithdrawalSerializer(many=True, required=False, allow_null=True)
    amount = serializers.ReadOnlyField()
    type = serializers.SlugRelatedField(slug_field="name", read_only=True)
    measurements = MeasurementSerializer(many=True, required=False, allow_null=True)

    class Meta:
        model = Well
        fields = "__all__"


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


class SimplePlateDimensionSerializer(serializers.ModelSerializer):
    class Meta:
        model = PlateDimension
        fields = ("id", "name", "rows", "cols")


class CompoundLibrarySerializer(serializers.ModelSerializer):
    class Meta:
        model = CompoundLibrary
        fields = "__all__"
        extra_kwargs = {
            "id": {
                "read_only": False,
                "required": False,
            },
        }


class PlateListSerializer(serializers.ModelSerializer):
    class Meta:
        model = Plate
        fields = "__all__"


class PlateDetailSerializer(serializers.ModelSerializer):
    class Meta:
        model = PlateDetail
        fields = "__all__"


class PlateSerializer(serializers.ModelSerializer):
    dimension = PlateDimensionSerializer(required=False, allow_null=True)
    details = serializers.SerializerMethodField()
    wells = serializers.SerializerMethodField()

    def get_details(self, plate: Plate):
        plate_details = PlateDetail.objects.get(pk=plate.id)
        return PlateDetailSerializer(plate_details).data

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

    # def to_representation(self, instance):
    #     representation = super().to_representation(instance)
    #     measurements_iso = instance.measurements

    #     for measurement in measurements_iso:
    #         measurement["measured_at"] = (
    #             measurement["measured_at"].isoformat().split("+")[0]
    #         )
    #     representation["measurements"] = measurements_iso
    #     return representation

    class Meta:
        model = Plate
        fields = "__all__"


class BarcodeSpecificationSerializer(serializers.ModelSerializer):
    class Meta:
        model = BarcodeSpecification
        fields = "__all__"


class ExperimentSerializer(serializers.ModelSerializer):
    plates = SimplePlateSerializer(many=True, required=False, allow_null=True)
    barcode_specifications = BarcodeSpecificationSerializer(
        many=True, required=False, allow_null=True
    )

    class Meta:
        model = Experiment
        fields = "__all__"


class SimpleExperimentSerializer(serializers.ModelSerializer):
    class Meta:
        model = Experiment
        fields = ("id", "name")


class SimplePlateTemplateSerializer(serializers.ModelSerializer):
    category = serializers.SlugRelatedField(slug_field="name", read_only=True)

    class Meta:
        model = PlateTemplate
        fields = ("id", "name", "category")


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

    class Meta:
        model = Project
        fields = "__all__"


class PlateMappingSerializer(UndefinedAffineModelSerializer):
    def save(self, **kwargs):
        try:
            if not self.validated_data["id"]:
                # Object is created: We apply the mapping to the source plate.
                mapping_file = self.validated_data["mapping_file"]
                from_column = self.validated_data["from_column"]
                to_column = self.validated_data["to_column"]
                amount_column = self.validated_data["amount_column"]
                delimiter = self.validated_data["delimiter"]
                quotechar = self.validated_data["quotechar"]
                sourcePlate = Plate.objects.get(pk=self.validated_data["source_plate"])
                targetPlate = Plate.objects.get(pk=self.validated_data["target_plate"])
                if mapping_file and from_column and to_column:
                    sourcePlate.map(
                        MappingList.from_csv(
                            mapping_file["name"],
                            from_column,
                            to_column,
                            amount_column,
                            delimiter,
                            quotechar,
                        ),
                        targetPlate,
                    )
                else:
                    sourcePlate.copy(targetPlate, self.validated_data["amount"])
            # TODO: What do we do on an update or delete??
            return super().save(**kwargs)
        except MappingError as ex:
            raise serializers.ValidationError({"detail": ex})

    class Meta:
        model = PlateMapping
        fields = "__all__"
