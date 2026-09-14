"""
Wells, their compounds, withdrawals and types.
"""

from rest_framework import serializers

from ..models import Well, WellCompound, WellWithdrawal
from .measurements import MeasurementSerializer
from .plates import SimplePlateSerializer


class WellCompoundSerializer(serializers.ModelSerializer):
    name = serializers.SlugRelatedField(
        slug_field="name", source="compound", read_only=True
    )

    structure = serializers.SlugRelatedField(
        slug_field="structure", source="compound", read_only=True
    )

    class Meta:
        model = WellCompound
        fields = "__all__"


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


class WellSerializer(serializers.ModelSerializer):
    hr_position = serializers.ReadOnlyField()
    current_info = serializers.ReadOnlyField()
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
