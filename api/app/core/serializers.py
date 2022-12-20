from rest_framework import serializers
from .models import Well, Plate, Measurement, PlateDimension
from compoundlib.models import CompoundLibrary, Compound


class CompoundSerializer(serializers.ModelSerializer):
    class Meta:
        model = Compound
        fields = ('name', 'identifier', 'structure')


class MeasurementSerializer(serializers.ModelSerializer):
    class Meta:
        model = Measurement
        fields = '__all__'


class WellSerializer(serializers.ModelSerializer):
    measurements = MeasurementSerializer(many=True)
    hr_position = serializers.ReadOnlyField()
    compounds = CompoundSerializer(many=True)

    class Meta:
        model = Well
        fields = '__all__'


class PlateDimensionSerializer(serializers.ModelSerializer):
    class Meta:
        model = PlateDimension
        fields = '__all__'


class SimplePlateDimensionSerializer(serializers.ModelSerializer):
    class Meta:
        model = PlateDimension
        fields = ('id', 'name', 'rows', 'cols')


class CompoundLibrarySerializer(serializers.ModelSerializer):
    class Meta:
        model = CompoundLibrary
        fields = '__all__'


class PlateSerializer(serializers.ModelSerializer):
    dimension = PlateDimensionSerializer()
    library = CompoundLibrarySerializer()
    wells = WellSerializer(many=True)

    # def get_wells(self, instance):
    #     wells = instance.wells.all().order_by('position')
    #     return WellSerializer(wells, many=True).data

    class Meta:
        model = Plate
        fields = '__all__'


class SimplePlateSerializer(serializers.ModelSerializer):
    dimension = serializers.SlugRelatedField(read_only=True, slug_field='name')

    class Meta:
        model = Plate
        fields = ('id', 'barcode', 'dimension')
